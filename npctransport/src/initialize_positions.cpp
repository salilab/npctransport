/**
 *  \file simulation_data.cpp
 *  \brief description.
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/initialize_positions.h>
#include <IMP/example/randomizing.h>
#include <IMP/example/optimizing.h>
#include <IMP/npctransport/particle_types.h>
#include <IMP/core/DistancePairScore.h>
#include <IMP/core/XYZR.h>
#include <IMP/base/Pointer.h>
#include <IMP/base/object_macros.h>
#include <IMP/core/RestraintsScoringFunction.h>
IMPNPCTRANSPORT_BEGIN_NAMESPACE

void initialize_positions(SimulationData *sd,
                          const ParticlePairsTemp &extra_links) {
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
  for (unsigned int i=0; i< extra_links.size(); ++i) {
    double d= core::XYZR(extra_links[i][0]).get_radius()
        +  core::XYZR(extra_links[i][1]).get_radius();
    IMP_NEW(core::HarmonicDistancePairScore, link,
            (d,sd->get_backbone_k(), "linker ps"));
    Pointer<Restraint> r
        = IMP::create_restraint(link.get(), extra_links[i]);
    rss.push_back(r);
  }
  IMP_NEW(core::RestraintsScoringFunction, sf, (rss));

  // Now optimize:
  int dump_interval = sd->get_rmf_dump_interval();
  sd->get_rmf_writer()->set_period(dump_interval * 100);// reduce output rate:
  example::optimize_balls(sd->get_diffusers()->get_particles(),
                          sf->get_restraints(),
                          sd->get_cpc()->get_pair_filters(),
                          OptimizerStates(1, sd->get_rmf_writer()),
                          PROGRESS);
  IMP_LOG(TERSE, "Initial energy is " << sd->get_m()->evaluate(false)
          << std::endl);
  sd->get_rmf_writer()->set_period(dump_interval);// restore output rate
  sd->get_rmf_writer()->update();

  // unpin previously unpinned fgs (= allow optimization)
  for (unsigned int i=0; i< previously_unpinned.size(); ++i) {
    previously_unpinned[i].set_coordinates_are_optimized(true);
  }
}


IMPNPCTRANSPORT_END_NAMESPACE
