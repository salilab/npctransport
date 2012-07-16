/**
 *  \file npctransport/rmf_links.h
 *  \brief description
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_RMF_LINKS_H
#define IMPNPCTRANSPORT_RMF_LINKS_H

#include "npctransport_config.h"
#include "SimulationData.h"
#include <IMP/rmf/atom_links.h>
#include <IMP/rmf/link_macros.h>
#include <IMP/WeakPointer.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

class IMPNPCTRANSPORTEXPORT HierarchyLoadLink:
  public rmf::HierarchyLoadLink {
 public:
  HierarchyLoadLink(RMF::FileConstHandle fh, Model *m);
};

class IMPNPCTRANSPORTEXPORT HierarchySaveLink:
  public rmf::HierarchySaveLink {
  WeakPointer<SimulationData> sd_;
  RMF::BallFactory bf_;
  RMF::ColoredFactory cf_;
  void do_add_recursive(Particle *root, Particle *p, RMF::NodeHandle cur);
  void do_save_node(Particle *p,
                    RMF::NodeHandle n,
                    unsigned int frame);
  std::pair<double, algebra::Vector3Ds> get_sites(core::ParticleType t) const ;
  // for testing without sd
  compatibility::map<core::ParticleType, std::pair<double,
      algebra::Vector3Ds > > sites_;
 public:
  HierarchySaveLink(RMF::FileHandle fh);
  // for testing
  void add_sites(core::ParticleType t,
                 double range,
                 algebra::Vector3Ds sites) {
    sites_[t]=std::make_pair(range, sites);
  }
};

// for testing
IMPNPCTRANSPORTEXPORT void add_sites(RMF::FileHandle fh,
                                     core::ParticleType t,
                                     double radius,
                                     algebra::Vector3Ds sites);

IMP_DECLARE_LINKERS(Hierarchy, hierarchy, hierarchies,
                    atom::Hierarchy,atom::Hierarchies,
                    atom::Hierarchy,atom::Hierarchies,
                    (RMF::FileHandle fh),
                    (RMF::FileConstHandle fh, Model *m));

IMPNPCTRANSPORT_END_NAMESPACE

#endif  /* IMPNPCTRANSPORT_RMF_LINKS_H */
