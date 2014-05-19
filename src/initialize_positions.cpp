/**
 *  \file simulation_data.cpp
 *  \brief description.
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/initialize_positions.h>
#include <IMP/npctransport/randomize_particles.h>
#include <IMP/npctransport/util.h>
#include <IMP/core/DistancePairScore.h>
#include <IMP/core/XYZR.h>
#include <IMP/atom/BrownianDynamics.h>
#include <IMP/Pointer.h>
#include <IMP/exception.h>
#include <IMP/object_macros.h>
#include <IMP/RAII.h>
#include <IMP/core/RestraintsScoringFunction.h>
#include <IMP/core/ConjugateGradients.h>
#include <IMP/core/MonteCarlo.h>
#include <IMP/core/SerialMover.h>
#include <IMP/core/BallMover.h>
#include <IMP/container/ClosePairContainer.h>
#include <IMP/container/ConsecutivePairContainer.h>
#include <IMP/core/rigid_bodies.h>
#include <IMP/core/SphereDistancePairScore.h>
#include <IMP/log_macros.h>
#include <IMP/log.h>
#include <IMP/container/ListSingletonContainer.h>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/scoped_ptr.hpp>
#include <IMP/scoped.h>
#include <IMP/PairPredicate.h>
#include <IMP/container/generic.h>
#include <IMP/internal/graph_utility.h>
#include <IMP/npctransport/SimulationData.h>
#include <IMP/raii_macros.h>
#include <IMP/log.h>
#include <IMP/flags.h>
//#include <IMP/example/optimizing.h>
#include <IMP/npctransport/internal/boost_main.h>
#include <cmath>

IMPNPCTRANSPORT_BEGIN_NAMESPACE
bool show_dependency_graph = false;
base::AddBoolFlag show_dep_adder("show_dependency_graph",
                                 "Show the dependency graph",
                                 &show_dependency_graph);

namespace {

//! adds a BallMover for each optimized particle in ps
//! with move size magnitude = particle_radius*2*scale_factor
//! and wraps it in a SerialMover in random order
/*
core::MonteCarloMover *create_serial_mover(const ParticlesTemp &ps,
                                           double scale_factor = 1.0) {
  core::MonteCarloMovers movers;
  for (unsigned int i = 0; i < ps.size(); ++i) {
    core::XYZR d = core::XYZR(ps[i]);
    double scale = d.get_radius() * 2.0 * scale_factor;
    if (d.get_coordinates_are_optimized()) {
      movers.push_back(new core::BallMover(ParticlesTemp(1, ps[i]), scale));
    }
  }
  std::random_shuffle(movers.begin(), movers.end());
  IMP_NEW(core::SerialMover, sm, (get_as<core::MonteCarloMoversTemp>(movers)));

  return sm.release();
  }*/

  /*
//! rescale the size of the move in all BallMovers inside everal serial mover
//! that is in the MonteCarlo obejct mc, by scale_factor
void rescale_move_size(core::MonteCarlo *mc, double scale_factor = 1.0) {
  using namespace core;
  IMP_ALWAYS_CHECK(scale_factor > 0, "move size scale factor must be positive",
                   IMP::ValueException);
  for (unsigned int i = 0; i < mc->get_number_of_movers(); ++i) {
    MonteCarloMover *mc_i = mc->get_mover(i);
    SerialMover *sm = dynamic_cast<SerialMover *>(mc_i);
    if (sm == nullptr) break;
    for (unsigned int j = 0; j < sm->get_movers().size(); ++j) {
      MonteCarloMover *sm_j = sm->get_movers()[j];
      BallMover *bm = dynamic_cast<BallMover *>(sm_j);
      if (bm == nullptr) break;
      bm->set_radius(bm->get_radius() * scale_factor);
    }
  }
  } */


/**
   create a Monte-Carlo object for partices ps, scored by isd,
   with a soft sphere pair score ssps, and move step size scaling move_scaling

   @param ps - the particles
   @param rs - restraints used for scoring the MC
   @param excluded - predicates for excluded volume restraints calc
   @param save - an optional saver to save optimization to file
   @param debug - if true, save is added as an optimizer state to the MC
   @param n_movers - number of movrs to apply in each MC step
 */
/*IMP::Pointer<core::MonteCarlo> create_mc(
    const ParticlesTemp &ps,
    const RestraintsTemp& rs,
    const PairPredicates excluded,
    rmf::SaveOptimizerState *save, bool debug,
    unsigned int n_movers = 1) {
  IMP_ALWAYS_CHECK(n_movers > 0, "number of MC movers must be >0",
                   ValueException);
  Model *m = ps[0]->get_model();
  IMP_NEW(core::MonteCarlo, mc, (m));
  //  mc->set_score_threshold(.1);
  if (debug && save) {
    mc->add_optimizer_state(save);
  }
  core::MonteCarloMovers movers;
  for (unsigned int i = 0; i < n_movers; ++i) {
    movers.push_back(create_serial_mover(ps));
  }
  mc->set_movers(movers);

  // Scoring
  IMP_NEW(core::IncrementalScoringFunction, isf, (ps, rs));
  IMP_NEW(core::SoftSpherePairScore, ssps, (10));
  // use special incremental support for the non-bonded part
  isf->add_close_pair_score(ssps, 0, ps, excluded);
  // TODO (what does this mean?)
  // we are special casing the nbl term for montecarlo, but using all for CG
  mc->set_incremental_scoring_function(isf);
  return mc;
}

*/


/** An RAII class for rescaling the x0 length of a LinearWellPairScore
    by some factor f (upon construction or using set()). Restores the
    original value upon destruction
*/
class LinearWellSetLengthRAII : public base::RAII {

  WeakPointer< LinearWellPairScore > ps_;
  double orig_;
  bool was_set_;

 public:
  IMP_RAII
    ( LinearWellSetLengthRAII, (LinearWellPairScore *ps, double f),
      { //Initialize
        was_set_ = false;
      },
      { // Set
        was_set_ = true;
        ps_ = ps;
        orig_ = ps_->get_rest_length_factor();
        ps_->set_rest_length_factor(orig_ * f);
      },
      { // Reset
        if(was_set_){
          ps_->set_rest_length_factor(orig_);
        }
      },
      { // Show });
      } );
};



/** An RAII class for temporarily changing the scoring function
    of an optimizer
*/
class OptimizerSetTemporaryScoringFunctionRAII: public base::RAII {

  PointerMember< IMP::ScoringFunction > orig_sf_;
  PointerMember< IMP::Optimizer > o_;
  bool was_set_;

 public:
  IMP_RAII
    ( OptimizerSetTemporaryScoringFunctionRAII,
      (Optimizer* o, ScoringFunctionAdaptor sfa),
      { //Initialize
        was_set_ = false;
      },
      { // Set
        was_set_ = true;
        o_ = o;
        orig_sf_ = o_->get_scoring_function();
        o_->set_scoring_function(sfa);
      },
      { // Reset
        if(was_set_){
          o_->set_scoring_function(orig_sf_);
        }
      },
      { // Show
      } );
};

/** An RAII class for temporarily changing the temperature
    of a BrownianDynamics object
*/
  class BDSetTemporaryTemperatureRAII : public base::RAII {

  double orig_temp_;
    base::PointerMember
    < IMP::atom::BrownianDynamics > bd_;
  bool was_set_;

 public:
  IMP_RAII
  ( BDSetTemporaryTemperatureRAII,
    (atom::BrownianDynamics* bd, double temp),
      { //Initialize
        was_set_ = false;
      },
      { // Set
        was_set_ = true;
        bd_ = bd;
        orig_temp_ = bd_->get_temperature();
        bd_->set_temperature(temp);
      },
      { // Reset
        if(was_set_){
          bd_->set_temperature(orig_temp_);
        }
      },
      { // Show
      } );
};

/** An RAII class for temporarily setting the optimization stsate of particles
*/
class TemporarySetOptimizationStateRAII: public base::RAII {
  bool orig_;
  WeakPointer<IMP::Model> m_;
  ParticleIndex pi_;
  bool was_set_;

 public:
  IMP_RAII
  ( TemporarySetOptimizationStateRAII,
    (ParticleAdaptor pa, bool is_optimized),
      { //Initialize
        was_set_ = false;
      },
      { // Set
        was_set_ = true;
        m_ = pa.get_model();
        pi_ = pa.get_particle_index();
        IMP_ALWAYS_CHECK(core::XYZ::get_is_setup( m_, pi_ ),
                         "p is not XYZ - can't set coordinates opt state",
                         IMP::ValueException);
        core::XYZ xyz( m_, pi_ );
        orig_ = xyz.get_coordinates_are_optimized();
        xyz.set_coordinates_are_optimized( is_optimized );
      },
      { // Reset
        if(was_set_){
          core::XYZ(m_, pi_).set_coordinates_are_optimized( orig_ );
        }
      },
      { // Show
      } );
};




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
                    LinearWellPairScores lwps,
                    base::LogLevel ll, bool debug,
                    double short_init_factor = 1.0) {
  // make sure that errors and log messages are marked as coming from this
  // function
  IMP_FUNCTION_LOG;

  base::SetLogState sls(ll);
  IMP_ALWAYS_CHECK(!ps.empty(), "No Particles passed.", ValueException);
  IMP_ALWAYS_CHECK(bd, "bd optimizer unspecified", ValueException);
  IMP_ALWAYS_CHECK(short_init_factor > 0 && short_init_factor <= 1.0,
                   "short init factor should be in range (0,1]",
                   ValueException);
  if(save)
    IMP_LOG(VERBOSE, "BEGIN o_b(): Saver has been called so far " <<
      save->get_number_of_updates() << std::endl);

  IMP_LOG(PROGRESS, "optimize_balls for n_particles="
          << ps.size() << std::endl);
  // shrink each of the particles, relax the configuration, repeat
  double bd_temperature_orig = bd->get_temperature();
  for (int i = 0; i < 11; ++i) {
    boost::scoped_array<boost::scoped_ptr<ScopedSetFloatAttribute> > attrs
      ( new boost::scoped_ptr<ScopedSetFloatAttribute>[ps.size()] );
    boost::ptr_vector<LinearWellSetLengthRAII> set_temporary_lengths;
    double factor = .1 * i;
    double length_factor = is_rest_length_scaling ? (.7 + .3 * factor) : 1.0;
    // rescale radii temporarily
    for (unsigned int j = 0; j < ps.size(); ++j) {
      attrs[j].reset(
          new ScopedSetFloatAttribute(ps[j], core::XYZR::get_radius_key(),
                                      core::XYZR(ps[j]).get_radius() * factor));
    }
    // rescale bond length temporarily
    for (unsigned int j = 0; j < lwps.size(); ++j) {
      set_temporary_lengths.push_back
        ( new LinearWellSetLengthRAII(lwps[j], length_factor) );
    }
    IMP_LOG(PROGRESS, "Optimizing with radii at " << factor << " of full"
            << " and length factor " << length_factor << std::endl
            << " energy before = "
            << bd->get_scoring_function()->evaluate(false) << std::endl);

    for (int k_simanneal = 0; k_simanneal < 5; ++k_simanneal)
      {
        double temperature =
          (1.5 - (i + k_simanneal / 5.0) / 11.0) * bd_temperature_orig;
        BDSetTemporaryTemperatureRAII
          bd_set_temporary_temperature(bd, temperature);
        bool done = false;
        IMP_OMP_PRAGMA(parallel num_threads(3)) {
          IMP_OMP_PRAGMA(single) {
            int n_bd_cycles =
              std::ceil(30 * (i / 2 + 2) * std::sqrt((double)ps.size()));
            int actual_n_bd_cycles =
              std::ceil(n_bd_cycles * short_init_factor) ;
            double e_bd = bd->optimize( actual_n_bd_cycles );
            IMP_LOG(PROGRESS, "Energy after bd is " << e_bd << " at " << i << ", "
                    << k_simanneal << std::endl);
            if (debug) {
              std::ostringstream oss;
              oss << "Init after " << i << " " << k_simanneal;
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


} // namespace {}



namespace {
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
    core::XYZ d(cur_chain->beads[0]);
    IMP_LOG(PROGRESS, "d # " << k << " = " << d << std::endl);
    IMP_LOG(PROGRESS, "is optimizable = " << d.get_coordinates_are_optimized()
            << std::endl);
  }
}
}

void initialize_positions(SimulationData *sd,
                          const RestraintsTemp &extra_restraints,
                          bool debug,
                          double short_init_factor) {
  IMP_FUNCTION_LOG;
  sd->set_was_used(true);
  IMP_ALWAYS_CHECK(short_init_factor > 0 && short_init_factor <= 1.0,
                   "short init factor should be in range (0,1]",
                   ValueException);
  randomize_particles(sd->get_diffusers()->get_particles(), sd->get_box());
  if (sd->get_rmf_sos_writer()) {
    sd->get_rmf_sos_writer()->update();
  }
  // pin first link of fgs, if not already pinned
  boost::ptr_vector<TemporarySetOptimizationStateRAII> chain_pins;
  atom::Hierarchies chains = sd->get_fg_chain_roots();
  for (unsigned int i = 0; i < chains.size(); ++i) {
    base::Pointer<FGChain> chain = get_fg_chain(chains[i]);
    chain_pins.push_back
      ( new TemporarySetOptimizationStateRAII
        (chain->beads[0], false) );
  }
  // core::XYZs previously_unpinned;
  // atom::Hierarchies chains = sd->get_fg_chains();
  // for (unsigned int i = 0; i < chains.size(); ++i) {

  //   if (core::XYZ(chains[i].get_child(0)).get_coordinates_are_optimized()) {
  //     previously_unpinned.push_back(core::XYZ(chains[i].get_child(0)));
  //     core::XYZ(chains[i].get_child(0)).set_coordinates_are_optimized(false);
  //   }
  // }

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
      Pointer<ScoringFunction> sf =
        sd->get_scoring()->get_custom_scoring_function
        ( extra_restraints,
          cur_particles,
          cur_optimizable_particles,
          false /* no non-bonded attr potentials yet */);
      OptimizerSetTemporaryScoringFunctionRAII
        set_temporary_scoring_function( sd->get_bd(), sf );
      optimize_balls(cur_particles,
                     false /*is_scale_rest_length*/,
                     sd->get_rmf_sos_writer(),
                     sd->get_bd(),
                     sd->get_scoring()->get_chain_scores(),
                     base::PROGRESS,
                     debug, short_init_factor);
    }
  {
    // optimize everything now
    ParticlesTemp particles =
      sd->get_diffusers()->get_particles(); // that should include obstacles
    ParticlesTemp optimizable_particles =
      get_optimizable_particles( particles );
    Pointer<ScoringFunction> sf =
      sd->get_scoring()->get_custom_scoring_function
      ( extra_restraints,
        particles,
        optimizable_particles ,
        false /* no non-bonded attr potential yet */ );
    OptimizerSetTemporaryScoringFunctionRAII
      set_temporary_scoring_function( sd->get_bd(), sf );
    optimize_balls(sd->get_diffusers()->get_particles(),
                   false /*scale rest length*/,
                   sd->get_rmf_sos_writer(),
                   sd->get_bd(),
                   sd->get_scoring()->get_chain_scores(),
                   base::PROGRESS,
                   debug, short_init_factor);
    IMP_LOG(VERBOSE, "Custom energy after initialization is "
              << sd->get_bd()->get_scoring_function()->evaluate(false)
              << std::endl);
  }
  IMP_LOG(VERBOSE, "Simulation energy after initialization is " <<
          sd->get_bd()->get_scoring_function()->evaluate(false)
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
