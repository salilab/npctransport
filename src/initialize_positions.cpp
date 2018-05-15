/**
 *  \file simulation_data.cpp
 *  \brief description.
 *
 *  Copyright 2007-2018 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/initialize_positions.h>
#include <IMP/npctransport/FGChain.h>
#include <IMP/npctransport/randomize_particles.h>
#include <IMP/npctransport/internal/initialize_positions_RAIIs.h>
#include <IMP/npctransport/internal/TAMDChain.h>
#include <IMP/npctransport/util.h>
#include <IMP/scoped.h>
#include <IMP/Restraint.h>
#include <IMP/ScoringFunction.h>
#include <IMP/atom/BrownianDynamics.h>
#include <IMP/atom/TAMDParticle.h>
#include <IMP/Pointer.h>
#include <IMP/exception.h>
#include <IMP/object_macros.h>
#include <IMP/RAII.h>
#include <IMP/core/RestraintsScoringFunction.h>
#include <IMP/core/ConjugateGradients.h>
#include <IMP/core/MonteCarlo.h>
#include <IMP/core/SerialMover.h>
#include <IMP/flags.h>
#include <IMP/log.h>
#include <IMP/log_macros.h>
//#include <IMP/internal/graph_utility.h>
//#include <IMP/core/RestraintsScoringFunction.h>
#include <IMP/core/BallMover.h>
#include <IMP/core/Hierarchy.h>
#include <IMP/core/rigid_bodies.h>
#include <IMP/core/XYZR.h>
#include <IMP/container/ListSingletonContainer.h>
#include <IMP/container/generic.h>
//#include <IMP/internal/graph_utility.h>
#include <IMP/npctransport/SimulationData.h>
#include <IMP/npctransport/internal/boost_main.h>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/scoped_ptr.hpp>
#include <cmath>
#include <vector>

IMPNPCTRANSPORT_BEGIN_NAMESPACE
// bool show_dependency_graph = false;
// AddBoolFlag show_dep_adder("show_dependency_graph",
//                                  "Show the dependency graph",
//                                  &show_dependency_graph);

namespace {

  // returns true if p is XYZ and TAMDParticle decorated
  bool is_xyz_tamd(IMP::Particle* p)
  {
    return IMP::atom::TAMDParticle::get_is_setup(p) &&
      IMP::core::XYZ::get_is_setup(p);
  }

  // update the coordinates of all TAMD particles from their
  // reference particles
  void update_tamd_particles_coords_from_refs(atom::Hierarchy root)
  {
    Particles ps;
    IMP::core::gather( root, is_xyz_tamd, std::back_inserter(ps));
    for(unsigned int i=0; i < ps.size(); i++) {
      IMP::atom::TAMDParticle t(ps[i]);
      t.update_coordinates_from_ref();
    }
  }

  /** Take a set of core::XYZR particles and relax them relative to a set of
      restraints. Excluded volume is handle separately, so don't include it
      in the passed list of restraints.
      (Note: optimize only particles that are flagged as 'optimized')

      @param ps particles over which to optimize
      @param is_rest_length_scaling if true - gradually increase the size of
              chains rest lengths during optimization
      @param short_init_factor a factor >0 and 1 for decreasing (or increasing if >1)
              the number of optimization cycles at each round
  */
  void optimize_balls(const ParticlesTemp &ps,
                      bool is_rest_length_scaling,
                      rmf::SaveOptimizerState *save,
                      atom::BrownianDynamics *bd,
                      FGChains chains,
                      LogLevel ll, bool debug,
                      double short_init_factor = 1.0) {
    // make sure that errors and log messages are marked as coming from this
    // function
    IMP_FUNCTION_LOG;

    SetLogState sls(ll);
    {
      IMP_ALWAYS_CHECK(!ps.empty(), "No Particles passed.",
                       IMP::ValueException);
      IMP_ALWAYS_CHECK(bd, "bd optimizer unspecified",
                       IMP::ValueException);
      IMP_ALWAYS_CHECK(short_init_factor > 0 ,
                       "short init factor should be positive",
                       IMP::ValueException);
    }
    if(save)
      IMP_LOG(VERBOSE, "BEGIN o_b(): Saver has been called so far " <<
              save->get_number_of_updates() << std::endl);

    std::set<ParticleIndex> ps_opt_set; // set of optimizable particles
    {
      ParticleIndexes ps_opt = get_particle_indexes
        ( get_optimizable_particles( ps ) );
      std::copy( ps_opt.begin(), ps_opt.end(),
                 std::inserter( ps_opt_set, ps_opt_set.end() ) );
    }
    IMP_LOG(PROGRESS, "optimize_balls for n_particles="
            << ps.size() <<
            " " << ps_opt_set.size() << " optimizables" << std::endl);
    // ramp temperature, particles radii, bond lengths, reset TAMD if needed
    // relax the configuration, repeat
    double bd_temperature_orig = bd->get_temperature();
    double bd_time_step_orig = bd->get_maximum_time_step();
    for ( double ramp_level = 0.001;
          ramp_level < 1.01 ;
          ramp_level = 1.4*ramp_level+.01 )
      {
        typedef boost::scoped_ptr<ScopedSetFloatAttribute> t_ptr_ScopedSetFloatAttribute;
        boost::scoped_array<t_ptr_ScopedSetFloatAttribute>
          tmp_set_radii( new t_ptr_ScopedSetFloatAttribute[ ps.size() ] );
        boost::ptr_vector< internal::FGChainScaleRestLengthRAII > tmp_scale_rest_length;
        boost::ptr_vector< internal::TAMDChainScaleKRAII > tmp_disable_tamd_k;
        double radius_factor_non_opt= std::pow(ramp_level,0.1); // radius scaling for obstacles / anchors
        double radius_factor_opt= ramp_level; // radius scaling for optimizables
        double rest_length_factor = (1.0/radius_factor_opt) * // radius_factor_opt*rest_length ~ 1.0 (so bond length is not affected by temporary scaling down of balls, only by the outcome of is_rest_length_scaling)
          (is_rest_length_scaling ? ( 0.5 + 0.5 * std::pow(ramp_level,0.75) ) : 1.0);
        // rescale particles radii temporarily
        for (unsigned int j = 0; j < ps.size(); ++j) {
          bool is_optimizable= ps_opt_set.count(ps[j]->get_index())==1;
          core::XYZR xyzr_j(ps[j]);
          double scaled_radius= xyzr_j.get_radius() * (is_optimizable ? radius_factor_opt : radius_factor_non_opt);
          tmp_set_radii[j].reset
            ( new ScopedSetFloatAttribute
              (ps[j], core::XYZR::get_radius_key(), scaled_radius) );
        }
        // rescale bond length + TAMD (if needed ) temporarily for all chains
        for (unsigned int j = 0; j < chains.size(); ++j) {
          //      std::cout<< "optimize balls CHAIN: " << chains[j] << std::endl;
          tmp_scale_rest_length.push_back
            ( new internal::FGChainScaleRestLengthRAII(chains[j],
                                                       rest_length_factor) );
          internal::TAMDChain* as_tamd =
            dynamic_cast< internal::TAMDChain* >(chains[j].get());
          if(as_tamd){
            tmp_disable_tamd_k.push_back
              ( new internal::TAMDChainScaleKRAII(as_tamd, 0.0) );
          }
        }
        IMP_LOG(PROGRESS, "Optimizing with radii at " << radius_factor_opt << " of full for optimizables"
                << " and length factor " << rest_length_factor << std::endl
                << " energy before = "
                << bd->get_scoring_function()->evaluate(false) << std::endl);
        // Simanneal:
        for (int k_simanneal = 0; k_simanneal < 5; ++k_simanneal)
          {
            double temperature =
              (2.5 - 2 * (10*ramp_level + k_simanneal / 4.0) / 11.0) * bd_temperature_orig;
            IMP::npctransport::internal::BDSetTemporaryTemperatureRAII
              bd_set_temporary_temperature(bd, temperature);
            double time_step =
              (22 - 21 * (10*ramp_level + k_simanneal / 4.0) / 11.0) * bd_time_step_orig;
            IMP::npctransport::internal::BDSetTemporaryTimeStepRAII
              bd_set_temporary_time_step(bd, time_step);
            bool done = false;
            IMP_OMP_PRAGMA(parallel num_threads(3)) {
              IMP_OMP_PRAGMA(single) {
                int n_bd_cycles =
                  std::ceil(30 * (5 * ramp_level + 2) * std::sqrt((double)ps_opt_set.size()));
                int actual_n_bd_cycles =
                  std::ceil(n_bd_cycles * short_init_factor) ;
                double e_bd = bd->optimize( actual_n_bd_cycles );
                IMP_LOG(PROGRESS,
                        "Energy after bd is " << e_bd <<
                        " at ramp level " << ramp_level << ", "
                        << k_simanneal << std::endl);
                if (debug) {
                  std::ostringstream oss;
                  oss << "Init after " << ramp_level << " " << k_simanneal;
                  if (save) {
                    bd->get_scoring_function()->evaluate(false);
                    save->update_always(oss.str());
                  }
                  IMP_LOG(VERBOSE, "updating RMF " << oss.str() << std::endl);
                }
              }
            }
            if (done) break;
          }  // for k_simanneal
      } // for i
    std::cout << "Energy after optimize_balls() using custom scoring function is " <<
      bd->get_scoring_function()->evaluate(false) << std::endl;
    if(save)
      IMP_LOG(VERBOSE, "END o_b(): Saver has been called so far " <<
              save->get_number_of_updates() << std::endl);
  } // optimize_balls()


  // print the first atoms of all the fgs in sd
  void print_fgs(IMP::npctransport::SimulationData &sd)
  {
    using namespace IMP;
    using atom::Hierarchy;
    using atom::Hierarchies;

    static int call_num = 0;
    IMP_LOG(PROGRESS, "INITIALIZE POSITIONS print_fgs() - Call # " << ++call_num
            << std::endl);

    Hierarchy root = sd.get_root();
    Hierarchies chains = sd.get_fg_chain_roots( );
    for (unsigned int k = 0; k < chains.size(); k++)
      {
        Pointer<FGChain> cur_chain = get_fg_chain(chains[k]);
        cur_chain->set_was_used(true);
        core::XYZ d(cur_chain->get_bead(0));
        IMP_LOG(PROGRESS, "d # " << k << " = " << d << std::endl);
        IMP_LOG(PROGRESS, "is optimizable = " << d.get_coordinates_are_optimized()
                << std::endl);
      }
  } // print_fgs()

  //! optimization of FGs at a nascent stage, which requires
  //! optimizing one at a time with various tweaks to e.g.
  //! obstacles radii
  //!
  //! Returns the FG particles processed (of type fg_type)
  ParticlesTemp initialize_positions_fg_at_a_time
  ( core::ParticleType fg_type,
    SimulationData* sd,
    ParticlesTemp obstacles,
    const RestraintsTemp &extra_restraints,
    bool debug,
    double short_init_factor)
  {
    ParticlesTemp beads= sd->get_beads();
    const bool is_rest_length_scaling(true);
    atom::Hierarchy cur_fg_root =  sd->get_root_of_type(fg_type);
    ParticlesTemp cur_fg_beads = atom::get_leaves( cur_fg_root );
    ParticlesTemp cur_particles = obstacles + cur_fg_beads;
    // pin all other particles:
    std::sort(beads.begin(), beads.end());
    std::sort(cur_particles.begin(), cur_particles.end());
    ParticlesTemp other_particles;
    std::set_difference(beads.begin(), beads.end(),
                        cur_particles.begin(), cur_particles.end(),
                        std::inserter(other_particles,
                                      other_particles.begin()) );
    boost::ptr_vector<
      IMP::npctransport::internal::TemporarySetOptimizationStateRAII> pin_others_particles;
    atom::Hierarchies chains = sd->get_fg_chain_roots();
    for (unsigned int i = 0; i < other_particles.size(); ++i) {
      pin_others_particles.push_back
        ( new IMP::npctransport::internal::TemporarySetOptimizationStateRAII
          (other_particles[i], false) );
    }
    // optimize
    std::cout <<  "Optimizing " <<  cur_particles.size()
              << " particles; adding chain type " << fg_type;
    ParticlesTemp cur_non_optimizable_beads =
      get_non_optimizable_particles( cur_particles);
    ParticlesTemp cur_optimizable_beads =
      get_optimizable_particles( cur_particles);
    IMP_LOG(VERBOSE, " ; " << cur_optimizable_beads.size()
            << " optimizable" << std::endl);
    if( cur_particles.size() * cur_optimizable_beads.size() == 0){
      IMP_LOG( WARNING, "No optimizable particles of type " << fg_type );
      return ParticlesTemp(); // noting to optimize if no (optimizable) particles
    }

    // switch local to sf till end of scope
    Pointer<ScoringFunction> sf =
      sd->get_scoring()->get_custom_scoring_function
      ( extra_restraints,
        get_particle_indexes(cur_non_optimizable_beads),
        cur_optimizable_beads,
        false /* no non-bonded attr potentials yet */);
    IMP::npctransport::internal::OptimizerSetTemporaryScoringFunctionRAII
      set_temporary_scoring_function( sd->get_bd(), sf );
    // inflate obstacles temporarily:
    boost::scoped_array<boost::scoped_ptr<ScopedSetFloatAttribute> >
      tmp_set_radii( new boost::scoped_ptr<ScopedSetFloatAttribute>[obstacles.size()] );
    for (unsigned int j = 0; j < obstacles.size(); ++j) {
      core::XYZR xyzr(obstacles[j]);
      double scaled_radius= xyzr.get_radius() * 2.0;
      tmp_set_radii[j].reset
        ( new ScopedSetFloatAttribute
          (obstacles[j], core::XYZR::get_radius_key(), scaled_radius) );
    }

    optimize_balls(cur_particles,
                   is_rest_length_scaling,
                   sd->get_rmf_sos_writer(),
                   sd->get_bd(),
                   sd->get_scoring()->get_fg_chains(),
                   PROGRESS,
                   debug, short_init_factor);
    update_tamd_particles_coords_from_refs( sd->get_root() );
    return cur_fg_beads;
  }

  //! optimization that is assumed to already involve most beads
  void initialize_positions_of_specific_beads
  ( SimulationData* sd,
    ParticlesTemp beads,
    const RestraintsTemp &extra_restraints,
    bool debug,
    double short_init_factor)
  {
    std::cout <<  "Optimizing " <<  beads.size()
              << " particles; " << std::endl;
    // Optimize everything now (including floats - kaps and inerts):
    ParticlesTemp non_optimizable_beads =
      npctransport::get_non_optimizable_particles( beads );
    ParticlesTemp optimizable_beads =
      npctransport::get_optimizable_particles( beads );
    IMP_LOG(VERBOSE, " ; " << optimizable_beads.size()
            << " optimizable" << std::endl);
    Pointer<ScoringFunction> sf =
      sd->get_scoring()->get_custom_scoring_function
      ( extra_restraints,
        get_particle_indexes(non_optimizable_beads),
        optimizable_beads,
        false /* no non-bonded attr potential yet */ );
    IMP::npctransport::internal::OptimizerSetTemporaryScoringFunctionRAII
      set_temporary_scoring_function( sd->get_bd(), sf );
    optimize_balls(sd->get_beads(),
                   false,
                   sd->get_rmf_sos_writer(),
                   sd->get_bd(),
                   sd->get_scoring()->get_fg_chains(),
                   PROGRESS,
                   debug,
                   short_init_factor);
    update_tamd_particles_coords_from_refs( sd->get_root() );
  }

} // namespace {}



void initialize_positions(SimulationData *sd,
                          const RestraintsTemp &extra_restraints,
                          bool debug,
                          double short_init_factor,
                          bool is_disable_randomize,
                          bool are_fgs_pre_initialized) {
  IMP_FUNCTION_LOG;
  sd->set_was_used(true);
  IMP_ALWAYS_CHECK(short_init_factor > 0,
                   "short init factor should be positive",
                   IMP::ValueException);
  if(!is_disable_randomize && !are_fgs_pre_initialized){
    randomize_particles(sd->get_beads(), sd->get_box()); // randomize before FGs are frozen
  }
  if (sd->get_rmf_sos_writer()) {
    sd->get_rmf_sos_writer()->update();
  }
  // pin first link of fgs, if not already pinned, or all fg beads if pre-initialized
  boost::ptr_vector<
    IMP::npctransport::internal::TemporarySetOptimizationStateRAII> chain_pins;
  atom::Hierarchies chains = sd->get_fg_chain_roots();
  for (unsigned int i = 0; i < chains.size(); ++i) {
    Pointer<FGChain> chain = get_fg_chain(chains[i]);
    unsigned int n_pinned=0;
    if(are_fgs_pre_initialized){
      n_pinned= chain->get_number_of_beads();
    }
    for(unsigned int j= 0; j<n_pinned; j++){
      chain_pins.push_back
        ( new IMP::npctransport::internal::TemporarySetOptimizationStateRAII
          (chain->get_bead(j), false) );
    }
  }
  if(!is_disable_randomize && are_fgs_pre_initialized){
    randomize_particles(sd->get_beads(), sd->get_box()); // randomize only now that FGs are frozen
  }


  int dump_interval = sd->get_rmf_dump_interval_frames();
  if (sd->get_rmf_sos_writer()) {
    // Now optimize:
    if (!debug) {
      sd->get_rmf_sos_writer()
          ->set_period(dump_interval * 100);  // reduce output rate:
      IMP_LOG(PROGRESS, "init: sos writer set dump interval period "
              << dump_interval * 100 << std::endl);
    } else {
      sd->get_rmf_sos_writer()->set_period(100);
      IMP_LOG(PROGRESS,  "init: sos writer set dump interval period "
              << 100 << std::endl);
    }
  }

  // optimize each FG separately using a temporary custom scoring function
  // (using RAII class OptimizerSetTemporaryScoringFunction)
  ParticlesTemp obstacles = sd->get_obstacle_particles();
  typedef boost::unordered_set<core::ParticleType> ParticleTypeSet;
  ParticleTypeSet const& fg_types = sd->get_fg_chain_types();
  ParticlesTemp processed_fg_beads;
  ParticlesTemp beads = sd->get_beads(); // that should include obstacles
  for(ParticleTypeSet::const_iterator
        it = fg_types.begin(); it != fg_types.end(); it++)
    {
      if(are_fgs_pre_initialized){
        break;
      }
      ParticlesTemp cur_fg_beads=
        initialize_positions_fg_at_a_time(*it, sd, obstacles,
                                          extra_restraints, debug,
                                          short_init_factor);
      processed_fg_beads += cur_fg_beads;
    }
  // Optimize with obstacles + all FGs:
  ParticlesTemp cur_particles = obstacles + processed_fg_beads;
  if(cur_particles.size()>0 && !are_fgs_pre_initialized) {
    initialize_positions_of_specific_beads(sd,
                                           cur_particles,
                                           extra_restraints,
                                           debug,
                                           short_init_factor*.5); // split with next step
  }
  // Final optimization of everything:
  initialize_positions_of_specific_beads(sd,
                                         beads,
                                         extra_restraints,
                                         debug,
                                         short_init_factor*.5);

  IMP_LOG(VERBOSE, "Custom energy after initialization is "
          << sd->get_bd()->get_scoring_function()->evaluate(true)
          << std::endl);

  IMP_LOG(VERBOSE, "Simulation energy after initialization is " <<
          sd->get_bd()->get_scoring_function()->evaluate(true)
          << std::endl);
  if (sd->get_rmf_sos_writer()) {
    sd->get_rmf_sos_writer()->set_period(dump_interval);  // restore output rate
    sd->get_bd()->get_scoring_function()->evaluate(false);
    IMP::RestraintsTemp rs(sd->get_scoring()->get_scoring_function_restraints(true));
    for(unsigned int i=0; i<rs.size(); i++){
      rs[i]->get_score();
    }
    sd->get_rmf_sos_writer()->update_always("done initializing");
    IMP_LOG(PROGRESS, "init: sos writer set dump interval period "
            << dump_interval << std::endl);
  }
  print_fgs(*sd);
}

IMPNPCTRANSPORT_END_NAMESPACE
