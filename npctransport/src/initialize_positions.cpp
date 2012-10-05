/**
 *  \file simulation_data.cpp
 *  \brief description.
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/initialize_positions.h>
#include <IMP/example/randomizing.h>
#include <IMP/npctransport/particle_types.h>
#include <IMP/core/DistancePairScore.h>
#include <IMP/core/XYZR.h>
#include <IMP/base/Pointer.h>
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
#include <IMP/npctransport/SimulationData.h>
#include <IMP/base/raii_macros.h>
#include <IMP/log.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

namespace {

core::Mover* create_serial_mover(const ParticlesTemp &ps) {
  core::Movers movers;
  for (unsigned int i=0; i< ps.size(); ++i) {
    double scale= core::XYZR(ps[i]).get_radius();
    movers.push_back(new core::BallMover(ParticlesTemp(1, ps[i]),
                                         scale*2));
  }
  IMP_NEW(core::SerialMover, sm, (get_as<core::MoversTemp>(movers)));
  return sm.release();
}

class SetLength: public base::RAII {
    base::WeakPointer<LinearWellPairScore> ps_;
    double orig_;
 public:
    IMP_RAII(SetLength, (LinearWellPairScore *ps, double f),
             {},
             {
               ps_=ps;
               orig_= ps_->get_x0();
               ps_->set_x0(orig_*f);
             },
             {
               ps_->set_x0(orig_);
             }, {});
  };


/** Take a set of core::XYZR particles and relax them relative to a set of
    restraints. Excluded volume is handle separately, so don't include it
in the passed list of restraints. */
void optimize_balls(const ParticlesTemp &ps,
                    const RestraintsTemp &rs,
                    const PairPredicates &excluded,
                    rmf::SaveOptimizerState *save,
                    Optimizer *local,
                    LinearWellPairScores dps,
                    base::LogLevel ll,
                    bool debug) {
  // make sure that errors and log messages are marked as coming from this
  // function
  IMP_FUNCTION_LOG;
  base::SetLogState sls(ll);
  IMP_USAGE_CHECK(!ps.empty(), "No Particles passed.");
  Model *m= ps[0]->get_model();
  //double scale = core::XYZR(ps[0]).get_radius();

  IMP_NEW(core::SoftSpherePairScore, ssps, (10));
  IMP_NEW(core::MonteCarlo, mc, (m));
  mc->set_score_threshold(.1);
  if (debug) {
    mc->add_optimizer_state(save);
  }
  IMP_NEW(core::IncrementalScoringFunction, isf, (ps, rs));
  {
    // set up MC
    mc->set_score_threshold(ps.size()*.1);
    mc->add_mover(create_serial_mover(ps));
    // we are special casing the nbl term for montecarlo, but using all for CG
    mc->set_incremental_scoring_function(isf);
    // use special incremental support for the non-bonded part
    isf->add_close_pair_score(ssps, 0, ps, excluded);
    // make pointer vector
  }

  IMP_LOG(PROGRESS, "Performing initial optimization" << std::endl);
  // shrink each of the particles, relax the configuration, repeat
  for (int i=0; i< 11; ++i) {
    boost::scoped_array<boost::scoped_ptr<ScopedSetFloatAttribute> > attrs
        (new boost::scoped_ptr<ScopedSetFloatAttribute>[ps.size()]);
    boost::ptr_vector<SetLength> lengths;
    double factor=.1*i;
    double length_factor=.3+.7*factor;
#pragma omp critical
    std::cout << "Optimizing with radii at " << factor << " of full"
              << " and length factor " << length_factor << std::endl;
    for (unsigned int j=0; j< ps.size(); ++j) {
      attrs[j].reset( new ScopedSetFloatAttribute(ps[j],
                                                  core::XYZR::get_radius_key(),
                                                  core::XYZR(ps[j])
                                                  .get_radius()*factor));
    }
    for (unsigned int j=0; j< dps.size(); ++j) {
      lengths.push_back(new SetLength(dps[j],
                                      length_factor));
    }
    // changed all radii
    isf->set_moved_particles(isf->get_movable_particles());
    for (int j=0; j< 5; ++j) {
      mc->set_kt(100.0/(3*j+1));
      double e= mc->optimize(ps.size()*(j+1)*500);
#pragma omp critical
      std::cout << "Energy is " << e << " at " << i << ", " << j << std::endl;
      if (debug) {
        std::ostringstream oss;
        oss << i << " " << j;
        save->update_always();
        save->set_frame_name(oss.str());
      }
      if (e < .000001) break;
    }
    double e=local->optimize(1000);
#pragma omp critical
    std::cout << "Energy after bd is " << e << std::endl;
  }
}
}

void initialize_positions(SimulationData *sd,
                          //                          const ParticlePairsTemp &extra_links,
                          const RestraintsTemp &extra_restraints,
                          bool debug) {
  example::randomize_particles(sd->get_diffusers()->get_particles(),
                               sd->get_box());
  sd->get_rmf_writer()->update();
  RestraintsTemp rss=sd->get_chain_restraints();
  if(sd->get_has_bounding_box())  rss.push_back(sd->get_box_restraint());
  if(sd->get_has_slab()) rss.push_back(sd->get_slab_restraint());
  // pin first link of fgs, if not already pinned
  core::XYZs previously_unpinned;
  atom::Hierarchies chains= get_fg_chains(sd->get_root());
  for (unsigned int i=0; i< chains.size(); ++i) {
    if (core::XYZ(chains[i].get_child(0)).get_coordinates_are_optimized()) {
      previously_unpinned.push_back(core::XYZ(chains[i].get_child(0)));
      core::XYZ(chains[i].get_child(0)).set_coordinates_are_optimized(false);
    }
  }
  rss += extra_restraints;
  // for (unsigned int i=0; i< extra_links.size(); ++i) {
  //   double d= core::XYZR(extra_links[i][0]).get_radius()
  //       +  core::XYZR(extra_links[i][1]).get_radius();
  //   IMP_NEW(core::HarmonicDistancePairScore, link,
  //           (d,sd->get_backbone_k(), "linker ps"));
  //   Pointer<Restraint> r
  //       = IMP::create_restraint(link.get(), extra_links[i]);
  //   rss.push_back(r);
  // }

  // Now optimize:
  int dump_interval = sd->get_rmf_dump_interval_frames();
  if (!debug) {
    sd->get_rmf_writer()->set_period(dump_interval * 100);// reduce output rate:
  } else {
    sd->get_rmf_writer()->set_period(100);
  }
  optimize_balls(sd->get_diffusers()->get_particles(),
                 rss,
                 sd->get_cpc()->get_pair_filters(),
                 sd->get_rmf_writer(),
                 sd->get_bd(),
                 sd->get_backbone_scores(),
                 PROGRESS, debug);
  IMP_NEW(core::RestraintsScoringFunction, rsf,
          (rss +RestraintsTemp(1, sd->get_predr()), "all restaints"));
  std::cout << "Initial energy is " << rsf->evaluate(false)
            << std::endl;
  sd->get_rmf_writer()->set_period(dump_interval);// restore output rate
  sd->get_rmf_writer()->update_always();
  sd->get_rmf_writer()->set_frame_name("done initializing");

  // unpin previously unpinned fgs (= allow optimization)
  for (unsigned int i=0; i< previously_unpinned.size(); ++i) {
    previously_unpinned[i].set_coordinates_are_optimized(true);
  }
}


IMPNPCTRANSPORT_END_NAMESPACE
