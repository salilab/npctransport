/**
 *  \file simulation_data.cpp
 *  \brief description.
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
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
#include <IMP/base/exception.h>
#include <IMP/base/flags.h>
#include <IMP/base/log.h>
#include <IMP/base/log_macros.h>
#include <IMP/base/object_macros.h>
#include <IMP/base/Pointer.h>
//#include <IMP/base/internal/graph_utility.h>
//#include <IMP/core/RestraintsScoringFunction.h>
#include <IMP/core/BallMover.h>
#include <IMP/core/Hierarchy.h>
#include <IMP/core/rigid_bodies.h>
#include <IMP/core/XYZR.h>
#include <IMP/container/ListSingletonContainer.h>
#include <IMP/container/generic.h>
#include <IMP/npctransport/SimulationData.h>
#include <IMP/npctransport/internal/boost_main.h>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/scoped_ptr.hpp>
#include <cmath>
#include <vector>

IMPNPCTRANSPORT_BEGIN_NAMESPACE
// bool show_dependency_graph = false;
// base::AddBoolFlag show_dep_adder("show_dependency_graph",
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
    @param short_init_factor a factor between >0 and 1 for decreasing
                             the number of optimization cycles at each
                             round
*/
void optimize_balls(const ParticlesTemp &ps,
                    bool is_rest_length_scaling,
                    rmf::SaveOptimizerState *save,
                    atom::BrownianDynamics *bd,
                    FGChains chains,
                    base::LogLevel ll, bool debug,
                    double short_init_factor = 1.0) {
  // make sure that errors and log messages are marked as coming from this
  // function
  IMP_FUNCTION_LOG;

  base::SetLogState sls(ll);
  {
    IMP_ALWAYS_CHECK(!ps.empty(), "No Particles passed.",
                     IMP::base::ValueException);
    IMP_ALWAYS_CHECK(bd, "bd optimizer unspecified",
                     IMP::base::ValueException);
    IMP_ALWAYS_CHECK(short_init_factor > 0 && short_init_factor <= 1.0,
                     "short init factor should be in range (0,1]",
                     IMP::base::ValueException);
  }
  if(save)
    IMP_LOG(VERBOSE, "BEGIN o_b(): Saver has been called so far " <<
      save->get_number_of_updates() << std::endl);

  IMP_LOG(PROGRESS, "optimize_balls for n_particles="
          << ps.size() << std::endl);

  // ramp temperature, particles radii, bond lengths, reset TAMD if needed
  // relax the configuration, repeat
  double bd_temperature_orig = bd->get_temperature();
  for (double ramp_level = 0; ramp_level <= 1.0 ; ramp_level += 0.1) {
    boost::scoped_array<boost::scoped_ptr<ScopedSetFloatAttribute> > tmp_set_radii
      ( new boost::scoped_ptr<ScopedSetFloatAttribute>[ps.size()] );
    boost::ptr_vector< internal::FGChainScaleRestLengthRAII >
      tmp_scale_rest_length;
    boost::ptr_vector< internal::TAMDChainScaleKRAII >
      tmp_disable_tamd_k;
    double radius_factor = ramp_level;
    double rest_length_factor =
      is_rest_length_scaling ? (.7 + .3 * ramp_level) : 1.0;
    // rescale particles radii temporarily
    for (unsigned int j = 0; j < ps.size(); ++j) {
      tmp_set_radii[j].reset(
          new ScopedSetFloatAttribute
          (ps[j], core::XYZR::get_radius_key(),
           core::XYZR(ps[j]).get_radius() * radius_factor));
    }
    // rescale bond length + TAMD (if needed ) temporarily for all chains
    for (unsigned int j = 0; j < chains.size(); ++j) {
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
    IMP_LOG(PROGRESS, "Optimizing with radii at " << radius_factor << " of full"
            << " and length factor " << rest_length_factor << std::endl
            << " energy before = "
            << bd->get_scoring_function()->evaluate(false) << std::endl);

    for (int k_simanneal = 0; k_simanneal < 5; ++k_simanneal)
      {
        double temperature =
          (1.5 - (10*ramp_level + k_simanneal / 5.0) / 11.0) * bd_temperature_orig;
        IMP::npctransport::internal::BDSetTemporaryTemperatureRAII
          bd_set_temporary_temperature(bd, temperature);
        bool done = false;
        IMP_OMP_PRAGMA(parallel num_threads(3)) {
          IMP_OMP_PRAGMA(single) {
            int n_bd_cycles =
              std::ceil(30 * (5 * ramp_level + 2) * std::sqrt((double)ps.size()));
            int actual_n_bd_cycles =
              std::ceil(n_bd_cycles * short_init_factor) ;
            double e_bd = bd->optimize( actual_n_bd_cycles );
            IMP_LOG(PROGRESS, "Energy after bd is " << e_bd <<
                    " at ramp level " << ramp_level << ", "
                    << k_simanneal << std::endl);
            if (debug) {
              std::ostringstream oss;
              oss << "Init after " << ramp_level << " " << k_simanneal;
              if (save) {
                save->update_always(oss.str());
              }
              IMP_LOG(VERBOSE, "updating RMF " << oss.str() << std::endl);
            }
          }
        }
        if (done) break;
      }  // for k_simanneal
  } // for i
  IMP_LOG(PROGRESS, "Energy after optimize_balls() is " <<
          bd->get_scoring_function()->evaluate(false) << std::endl);
  if(save)
    IMP_LOG(VERBOSE, "END o_b(): Saver has been called so far " <<
      save->get_number_of_updates() << std::endl);
}


// print the first atoms of all the fgs in sd
void print_fgs(IMP::npctransport::SimulationData &sd) {
  using namespace IMP;
  using atom::Hierarchy;
  using atom::Hierarchies;

  static int call_num = 0;
  IMP_LOG(PROGRESS, "INITIALIZE POSITIONS print_fgs() - Call # " << ++call_num
          << std::endl);

  Hierarchy root = sd.get_root();
  Hierarchies chains = sd.get_fg_chain_roots( );
  for (unsigned int k = 0; k < chains.size(); k++) {
    base::Pointer<FGChain> cur_chain = get_fg_chain(chains[k]);
    core::XYZ d(cur_chain->get_bead(0));
    IMP_LOG(PROGRESS, "d # " << k << " = " << d << std::endl);
    IMP_LOG(PROGRESS, "is optimizable = " << d.get_coordinates_are_optimized()
            << std::endl);
  }
}

} // namespace {}


void initialize_positions(SimulationData *sd,
                          const RestraintsTemp &extra_restraints,
                          bool debug,
                          double short_init_factor) {
  IMP_FUNCTION_LOG;
  sd->set_was_used(true);
  IMP_ALWAYS_CHECK(short_init_factor > 0 && short_init_factor <= 1.0,
                   "short init factor should be in range (0,1]",
                   IMP::base::ValueException);
  randomize_particles(sd->get_diffusers()->get_particles(), sd->get_box());
  if (sd->get_rmf_sos_writer()) {
    sd->get_rmf_sos_writer()->update();
  }
  // pin first link of fgs, if not already pinned
  boost::ptr_vector<
    IMP::npctransport::internal::TemporarySetOptimizationStateRAII> chain_pins;
  atom::Hierarchies chains = sd->get_fg_chain_roots();
  for (unsigned int i = 0; i < chains.size(); ++i) {
    base::Pointer<FGChain> chain = get_fg_chain(chains[i]);
    chain_pins.push_back
      ( new IMP::npctransport::internal::TemporarySetOptimizationStateRAII
        (chain->get_bead(0), false) );
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
  ParticleTypeSet const& types = sd->get_fg_types();
  ParticlesTemp cur_particles = obstacles;
  for(ParticleTypeSet::const_iterator
        it = types.begin(); it != types.end(); it++)
    {
      atom::Hierarchy cur_fg_root =  sd->get_root_of_type(*it);
      ParticlesTemp cur_fg_beads = atom::get_leaves( cur_fg_root );
      cur_particles += cur_fg_beads;
      IMP_LOG(VERBOSE, "Optimizing " <<  cur_particles.size()
               << " particles of type " << *it);
      ParticlesTemp cur_optimizable_particles =
        get_optimizable_particles( cur_particles);
      IMP_LOG(VERBOSE, " ; " << cur_optimizable_particles.size()
             << " optimizable" << std::endl);
      if( cur_particles.size() * cur_optimizable_particles.size() == 0){
        IMP_LOG( WARNING, "No optimizable particles of type " << *it );
        continue; // noting to optimize if no (optimizable) particles
      }
      // switch local to sf till end of scope
      base::Pointer<ScoringFunction> sf =
        sd->get_scoring()->get_custom_scoring_function
        ( extra_restraints,
          cur_particles,
          cur_optimizable_particles,
          false /* no non-bonded attr potentials yet */);
      IMP::npctransport::internal::OptimizerSetTemporaryScoringFunctionRAII
        set_temporary_scoring_function( sd->get_bd(), sf );
      optimize_balls(cur_particles,
                     false /*is_scale_rest_length*/,
                     sd->get_rmf_sos_writer(),
                     sd->get_bd(),
                     sd->get_scoring()->get_fg_chains(),
                     base::PROGRESS,
                     debug, short_init_factor);
      update_tamd_particles_coords_from_refs( sd->get_root() );
    }
  {
    // optimize everything now
    ParticlesTemp particles =
      sd->get_diffusers()->get_particles(); // that should include obstacles
    ParticlesTemp optimizable_particles =
      get_optimizable_particles( particles );
    base::Pointer<ScoringFunction> sf =
      sd->get_scoring()->get_custom_scoring_function
      ( extra_restraints,
        particles,
        optimizable_particles ,
        false /* no non-bonded attr potential yet */ );
    IMP::npctransport::internal::OptimizerSetTemporaryScoringFunctionRAII
      set_temporary_scoring_function( sd->get_bd(), sf );
    optimize_balls(sd->get_diffusers()->get_particles(),
                   false /*scale rest length*/,
                   sd->get_rmf_sos_writer(),
                   sd->get_bd(),
                   sd->get_scoring()->get_fg_chains(),
                   base::PROGRESS,
                   debug, short_init_factor);
    update_tamd_particles_coords_from_refs( sd->get_root() );
    IMP_LOG(VERBOSE, "Custom energy after initialization is "
              << sd->get_bd()->get_scoring_function()->evaluate(true)
              << std::endl);
  }
  IMP_LOG(VERBOSE, "Simulation energy after initialization is " <<
          sd->get_bd()->get_scoring_function()->evaluate(true)
          << std::endl);
  print_fgs(*sd);
  // IMP_NEW(core::RestraintsScoringFunction, rsf,
  //         (rss + RestraintsTemp
  //          ( 1, sd->get_scoring()->get_predicates_pair_restraint() ),
  //          "all restaints"));
  // IMP_LOG(WARNING, "Initial energy is " << rsf->evaluate(false) << std::endl);
  if (sd->get_rmf_sos_writer()) {
    sd->get_rmf_sos_writer()->set_period(dump_interval);  // restore output rate
    sd->get_rmf_sos_writer()->update_always("done initializing");
    IMP_LOG(base::PROGRESS, "init: sos writer set dump interval period "
            << dump_interval << std::endl);
  }
  print_fgs(*sd);
}

IMPNPCTRANSPORT_END_NAMESPACE
