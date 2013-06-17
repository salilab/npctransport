/**
 *  \file simulation_data.cpp
 *  \brief description.
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/initialize_positions.h>
#include <IMP/npctransport/randomize_particles.h>
#include <IMP/npctransport/particle_types.h>
#include <IMP/core/DistancePairScore.h>
#include <IMP/core/XYZR.h>
#include <IMP/atom/BrownianDynamics.h>
#include <IMP/base/Pointer.h>
#include <IMP/base/exception.h>
#include <IMP/base/object_macros.h>
#include <IMP/core/RestraintsScoringFunction.h>
#include <IMP/core/ConjugateGradients.h>
#include <IMP/core/MonteCarlo.h>
#include <IMP/core/SerialMover.h>
#include <IMP/core/BallMover.h>
#include <IMP/container/ClosePairContainer.h>
#include <IMP/core/rigid_bodies.h>
#include <IMP/core/SphereDistancePairScore.h>
#include <IMP/base/log_macros.h>
#include <IMP/container/ListSingletonContainer.h>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/scoped_ptr.hpp>
#include <IMP/scoped.h>
#include <IMP/PairPredicate.h>
#include <IMP/container/generic.h>
#include <IMP/base/internal/graph_utility.h>
#include <IMP/npctransport/SimulationData.h>
#include <IMP/base/raii_macros.h>
#include <IMP/base/log.h>
#include <IMP/base/flags.h>
#include <IMP/example/optimizing.h>
#ifdef IMP_NPC_GOOGLE
#include <IMP/npctransport/internal/google_main.h>
#else
#include <IMP/npctransport/internal/boost_main.h>
#endif
IMPNPCTRANSPORT_BEGIN_NAMESPACE
bool short_initialize = false;
base::AddBoolFlag short_adder("short_initialize",
                              "Run an abbreviated version of initialize",
                              &short_initialize);
bool show_dependency_graph = false;
base::AddBoolFlag show_dep_adder("show_dependency_graph",
                                 "Show the dependency graph",
                                 &show_dependency_graph);

namespace {

//! adds a BallMover for each optimized particle in ps
//! with move size magnitude = particle_radius*2*scale_factor
//! and wraps it in a SerialMover in random order
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
}

//! rescale the size of the move in all BallMovers inside everal serial mover
//! that is in the MonteCarlo obejct mc, by scale_factor
void rescale_move_size(core::MonteCarlo *mc, double scale_factor = 1.0) {
  using namespace core;
  IMP_ALWAYS_CHECK(scale_factor > 0, "move size scale factor must be positive",
                   IMP::base::ValueException);
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
}

class SetLength : public base::RAII {
  base::WeakPointer<LinearWellPairScore> ps_;
  double orig_;

 public:
  IMP_RAII(SetLength, (LinearWellPairScore *ps, double f), {}, {
    ps_ = ps;
    orig_ = ps_->get_x0();
    ps_->set_x0(orig_ * f);
  },
           {
    ps_->set_x0(orig_);
  },
           {});
};

/**
   create a Monte-Carlo object for partices ps, scored by isd,
   with a soft sphere pair score ssps, and move step size scaling move_scaling

   @param ps - the particles
   @param excluded - predicates for excluded volume restraints calc
   @param save - an optional saver to save optimization to file
   @param debug - if true, save is added as an optimizer state to the MC
   @param n_movers - number of movrs to apply in each MC step
 */
IMP::base::Pointer<core::MonteCarlo> create_mc(
    const ParticlesTemp &ps,
    IMP::base::Pointer<core::IncrementalScoringFunction> isf,
    IMP::base::Pointer<core::SoftSpherePairScore> ssps,
    const PairPredicates &excluded, rmf::SaveOptimizerState *save, bool debug,
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
  // we are special casing the nbl term for montecarlo, but using all for CG
  mc->set_incremental_scoring_function(isf);
  // use special incremental support for the non-bonded part
  isf->add_close_pair_score(ssps, 0, ps, excluded);
  return mc;
}

/** Take a set of core::XYZR particles and relax them relative to a set of
    restraints. Excluded volume is handle separately, so don't include it
    in the passed list of restraints.
    (Note: optimize only particles that are flagged as 'optimized')
*/
void optimize_balls(const ParticlesTemp &ps, const RestraintsTemp &rs,
                    const PairPredicates &excluded,
                    rmf::SaveOptimizerState *save,
                    atom::BrownianDynamics *local, LinearWellPairScores dps,
                    base::LogLevel ll, bool debug) {
  // make sure that errors and log messages are marked as coming from this
  // function
  IMP_FUNCTION_LOG;

  base::SetLogState sls(ll);
  IMP_ALWAYS_CHECK(!ps.empty(), "No Particles passed.", ValueException);
  IMP_ALWAYS_CHECK(local, "local optimizer unspecified", ValueException);
  Model *m = ps[0]->get_model();
  // double scale = core::XYZR(ps[0]).get_radius();

  IMP_NEW(core::IncrementalScoringFunction, isf, (ps, rs));
  IMP_NEW(core::SoftSpherePairScore, ssps, (10));
  unsigned int n_movers = 1 + ps.size() / 600;  // movers applied in one MC step
  IMP::base::Pointer<core::MonteCarlo> mc =
      create_mc(ps, isf, ssps, excluded, save, debug, n_movers);

  if (show_dependency_graph) {
    DependencyGraph dg = get_dependency_graph(ps[0]->get_model());
    std::ofstream dgf("dependency_graph.dot");
    show_as_graphviz(dg, dgf);
  }

  IMP_LOG(PROGRESS, "Performing initial optimization" << std::endl);
  // shrink each of the particles, relax the configuration, repeat
  double bd_temperature_orig = local->get_temperature();
  for (int i = 0; i < 11; ++i) {
    boost::scoped_array<boost::scoped_ptr<ScopedSetFloatAttribute> > attrs(
        new boost::scoped_ptr<ScopedSetFloatAttribute>[ps.size()]);
    boost::ptr_vector<SetLength> lengths;
    double factor = .1 * i;
    double length_factor = .7 + .3 * factor;
    // rescale radii
    for (unsigned int j = 0; j < ps.size(); ++j) {
      attrs[j].reset(
          new ScopedSetFloatAttribute(ps[j], core::XYZR::get_radius_key(),
                                      core::XYZR(ps[j]).get_radius() * factor));
    }
    // rescale bond length
    for (unsigned int j = 0; j < dps.size(); ++j) {
      lengths.push_back(new SetLength(dps[j], length_factor));
    }
    std::cout << "Optimizing with radii at " << factor << " of full"
              << " and length factor " << length_factor << std::endl;
    isf->set_moved_particles(isf->get_movable_particles());
    double desired_accept_rate =
        0.2;  // acceptance rate for which to tune move size (per mover)
    for (int j = 0; j < 5; ++j) {
      local->set_temperature((1.5 - (i + j / 5.0) / 11.0) *
                             bd_temperature_orig);
      double kt = 100.0 / (3 * j + 1);
      bool done = false;
      IMP_OMP_PRAGMA(parallel num_threads(3)) {
        int timer = 1000;
        IMP_OMP_PRAGMA(single) {
          int n_bd = 1000 * (i + 1);
          int n_bd_inner = 1000;
          int n_external = n_bd / n_bd_inner;
          int mc_opt_factor = 1 * (j + 1);  // 500
          int n_mc_all = ps.size() * mc_opt_factor / n_movers;
          int n_mc_inner = n_mc_all / n_external;
          for (int k = 0; k < n_external; k++) {
            // mc->set_kt(kt);
            // std::cout << "n_mc_inner = " << n_mc_inner<< std::endl;
            // double e_mc = mc->optimize(n_mc_inner);
            // unsigned int n_accepted =
            //   mc->get_mover(0)->get_number_of_accepted();
            // unsigned int n_proposed =
            //   mc->get_mover(0)->get_number_of_proposed();
            // double accept_rate = n_accepted / (n_proposed + 0.0) ;
            // std::cout << "Energy is " << e_mc << " at " << i <<
            //   ", " << j << "," << k << std::endl;
            // std::cout << "For each of " << n_movers << " movers, accepted "
            //           << n_accepted << " of " << n_proposed
            //           << " = " << accept_rate  << std::endl;
            // // scale move size to get 'desired_accept_rate' acceptance rate
            // double move_rescale;
            // if(accept_rate > desired_accept_rate) {
            //   move_rescale = accept_rate / desired_accept_rate;
            // } else {
            //   move_rescale =
            //     (accept_rate + 0.01*desired_accept_rate)
            //     / (accept_rate + 0.1 *desired_accept_rate) ;
            // }
            // std::cout << "Rescaling move size and kt by " << move_rescale <<
            // std::endl;
            // rescale_move_size( mc, move_rescale);
            // kt /= move_rescale;
            double e_bd = local->optimize(n_bd_inner);
            std::cout << "Energy after bd is " << e_bd << " at " << i << ", "
                      << j << "," << k << std::endl;
            //            mc->get_mover(0)->reset_statistics();
            if (--timer == 0) {/*exit(0);*/
            }
          }  // for k

          //          double e= mc->optimize(!short_initialize ?
          // ps.size()*(j+1)*500: 1000);
          // std::cout << "Energy is " << e << " at " << i << ", " << j <<
          // std::endl;
          if (debug) {
            std::ostringstream oss;
            oss << "Init after " << i << " " << j;
            if (save) {
              save->update_always(oss.str());
            }
          }
          // if (e < .000001) done=true; // TODO: replace stopping condition
        }
      }
      if (short_initialize) return;
      if (done) break;
    }  // for j
    IMP_OMP_PRAGMA(parallel num_threads(3)) {
      IMP_OMP_PRAGMA(single) {
        double e = local->optimize(1000);
        std::cout << "Energy after bd is " << e << std::endl;
      }
    }
  }
  local->set_temperature(bd_temperature_orig);
}
}

namespace {
// print the first atoms of all the fgs in sd
void print_fgs(IMP::npctransport::SimulationData &sd) {
  using namespace IMP;
  using atom::Hierarchy;
  using atom::Hierarchies;

  static int call_num = 0;
  std::cout << "INITIALIZE POSITIONS print_fgs() - Call # " << ++call_num
            << std::endl;

  Hierarchy root = sd.get_root();
  Hierarchies chains = IMP::npctransport::get_fg_chains(root);
  for (unsigned int k = 0; k < chains.size(); k++) {
    Hierarchy cur_chain(chains[k]);
    core::XYZ d(cur_chain.get_child(0));
    std::cout << "d # " << k << " = " << d << std::endl;
    std::cout << "is optimizable = " << d.get_coordinates_are_optimized()
              << std::endl;
  }
}
}

void initialize_positions(SimulationData *sd,
                          //                          const ParticlePairsTemp
                          // &extra_links,
                          const RestraintsTemp &extra_restraints,
                          bool debug) {
  sd->set_was_used(true);
  print_fgs(*sd);
  randomize_particles(sd->get_diffusers()->get_particles(), sd->get_box());
  print_fgs(*sd);
  if (sd->get_rmf_sos_writer()) {
    sd->get_rmf_sos_writer()->update();
  }
  RestraintsTemp rss = sd->get_chain_restraints();
  if (sd->get_has_bounding_box()) rss.push_back(sd->get_box_restraint());
  if (sd->get_has_slab()) rss.push_back(sd->get_slab_restraint());
  // pin first link of fgs, if not already pinned
  core::XYZs previously_unpinned;
  atom::Hierarchies chains = get_fg_chains(sd->get_root());
  for (unsigned int i = 0; i < chains.size(); ++i) {
    if (core::XYZ(chains[i].get_child(0)).get_coordinates_are_optimized()) {
      previously_unpinned.push_back(core::XYZ(chains[i].get_child(0)));
      core::XYZ(chains[i].get_child(0)).set_coordinates_are_optimized(false);
    }
  }
  print_fgs(*sd);
  rss += extra_restraints;
  // for (unsigned int i=0; i< extra_links.size(); ++i) {
  //   double d= core::XYZR(extra_links[i][0]).get_radius()
  //       +  core::XYZR(extra_links[i][1]).get_radius();
  //   IMP_NEW(core::HarmonicDistancePairScore, link,
  //           (d,sd->get_backbone_k(), "linker ps"));
  //   base::Pointer<Restraint> r
  //       = IMP::create_restraint(link.get(), extra_links[i]);
  //   rss.push_back(r);
  // }

  int dump_interval = sd->get_rmf_dump_interval_frames();
  if (sd->get_rmf_sos_writer()) {
    // Now optimize:
    if (!debug) {
      sd->get_rmf_sos_writer()
          ->set_period(dump_interval * 100);  // reduce output rate:
    } else {
      sd->get_rmf_sos_writer()->set_period(100);
    }
  }
  optimize_balls(sd->get_diffusers()->get_particles(), rss,
                 sd->get_cpc()->get_pair_filters(), sd->get_rmf_sos_writer(),
                 sd->get_bd(), sd->get_backbone_scores(), base::PROGRESS,
                 debug);
  print_fgs(*sd);
  IMP_NEW(core::RestraintsScoringFunction, rsf,
          (rss + RestraintsTemp(1, sd->get_predr()), "all restaints"));
  std::cout << "Initial energy is " << rsf->evaluate(false) << std::endl;
  if (sd->get_rmf_sos_writer()) {
    sd->get_rmf_sos_writer()->set_period(dump_interval);  // restore output rate
    sd->get_rmf_sos_writer()->update_always("done initializing");
  }
  // unpin previously unpinned fgs (= allow optimization)
  for (unsigned int i = 0; i < previously_unpinned.size(); ++i) {
    previously_unpinned[i].set_coordinates_are_optimized(true);
  }
  print_fgs(*sd);
}

IMPNPCTRANSPORT_END_NAMESPACE
