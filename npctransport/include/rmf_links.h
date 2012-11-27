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

class IMPNPCTRANSPORTEXPORT HierarchyWithSitesLoadLink:
  public rmf::HierarchyLoadLink {
    WeakPointer<SimulationData> sd_;
  RMF::BallConstFactory bf_;
  IMP_PROTECTED_METHOD(virtual void, do_add_link_recursive,
      (Particle *root,
                             Particle *o, RMF::NodeConstHandle node),
    IMP_OVERRIDE,);


 public:
  HierarchyWithSitesLoadLink(RMF::FileConstHandle fh, Model *m);
};

class IMPNPCTRANSPORTEXPORT HierarchyWithSitesSaveLink:
  public rmf::HierarchySaveLink {
  WeakPointer<SimulationData> sd_;
  RMF::BallFactory bf_;
  RMF::ColoredFactory cf_;
  IMP_PROTECTED_METHOD(virtual void, do_add_recursive,
              (Particle *root, Particle *p,
                             RMF::NodeHandle cur), IMP_OVERRIDE,);
  std::pair<double, algebra::Vector3Ds> get_sites(core::ParticleType t) const ;
  // for testing without sd
  compatibility::map<core::ParticleType, std::pair<double,
      algebra::Vector3Ds > > sites_;
 public:
  HierarchyWithSitesSaveLink(RMF::FileHandle fh);
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

// note that the corresponding define macro in the .cpp file implicitly
// uses the HierarchyWithSitesLoadLink and HierarchyWithSitesSaveLink
// classes
/**
   Functions for adding a hierarchy to an RMF file or linking an RMF
   file to an existing hierarchy, including support for particles with sites.

   These functions are practically identical to the add / link
   methods in modules/rmf/include/hierarchy_io.h, such as
   IMP::rmf::link_hierarchies(), except here NPC particle sites
   support is included.
 */
IMP_DECLARE_LINKERS(HierarchyWithSites,
                    hierarchy_with_sites, hierarchies_with_sites,
                    atom::Hierarchy,atom::Hierarchies,
                    atom::Hierarchy,atom::Hierarchies,
                    (RMF::FileConstHandle fh, Model *m),
                    See IMP::rmf::link_hierarchies() for more details.
                    The only difference is the addition of particle sites
                    support.
                    );

IMPNPCTRANSPORT_END_NAMESPACE

#endif  /* IMPNPCTRANSPORT_RMF_LINKS_H */
