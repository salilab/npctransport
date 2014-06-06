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
#include <boost/unordered_map.hpp>
#include <IMP/WeakPointer.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

class IMPNPCTRANSPORTEXPORT HierarchyWithSitesLoadLink
    : public rmf::HierarchyLoadLink {
  WeakPointer<SimulationData> sd_;
  RMF::decorator::BallConstFactory bf_;
  // Transporting decorator related keys:
  RMF::IntKey is_last_entry_from_top_key_;
  RMF::IntKey n_entries_bottom_key_;
  RMF::IntKey n_entries_top_key_;
  // optimization data:
  RMF::IntKey coordinates_are_optimized_key_;
  boost::unordered_map<ParticleIndex, ParticleIndexes> particles_;

 protected:
  /** link the specified rmf node to particle cur. In addition, load
      the list of sites associated with the rmf node together with
      their coordinates (children of type RMF::Geometry decorated with
      RMF::BallConst) and associate them with the particle type of cur
      in the SimulationData of cur.
  */
  virtual void do_link_particle(kernel::Model *m, kernel::ParticleIndex root,
                                kernel::ParticleIndex cur,
                                RMF::NodeConstHandle node) IMP_OVERRIDE;

  virtual void do_setup_particle(Model* m, ParticleIndex root,
                                 ParticleIndex cur,
                                 RMF::NodeConstHandle node) IMP_OVERRIDE {
    IMP_NOT_IMPLEMENTED;
    IMP_UNUSED(m);     IMP_UNUSED(root);
    IMP_UNUSED(cur);    IMP_UNUSED(node);
  }

  /** load the values of a hierarchy. also
      loads dynamic Transporting decorator transport directionality
      information if needed (which is used in
      ParticleTransportStatisticsOptimizerState)
  */
  virtual void do_load_hierarchy(RMF::NodeConstHandle root_node,
                                 Model *m,
                                 ParticleIndex pi) IMP_OVERRIDE;

 public:
  HierarchyWithSitesLoadLink(RMF::FileConstHandle fh);
  static const char *get_name() {return "npctransport load";}
};

class IMPNPCTRANSPORTEXPORT HierarchyWithSitesSaveLink
    : public rmf::HierarchySaveLink {
  WeakPointer<SimulationData> sd_;
  RMF::decorator::BallFactory bf_;
  RMF::decorator::ColoredFactory cf_;
  // Transporting decorator related keys:
  RMF::IntKey is_last_entry_from_top_key_;
  RMF::IntKey n_entries_bottom_key_;
  RMF::IntKey n_entries_top_key_;
  // optimization data:
  RMF::IntKey coordinates_are_optimized_key_;

  boost::unordered_map<ParticleIndex, ParticleIndexes> particles_;

  std::pair<double, algebra::Vector3Ds> get_sites(core::ParticleType t) const;
  // for testing without sd
  boost::unordered_map<core::ParticleType, std::pair<double, algebra::Vector3Ds> > sites_;

 protected:
  virtual void do_setup_node(Model *m, ParticleIndex root,
                           ParticleIndex cur,
                           RMF::NodeHandle cur_node) IMP_OVERRIDE;

  /** save the values of a hierarchy
  */
  virtual void do_save_hierarchy(Model *m, ParticleIndex root,
                                 RMF::NodeHandle root_node) IMP_OVERRIDE;

 public:
  HierarchyWithSitesSaveLink(RMF::FileHandle fh);
  // for testing
  void add_sites(core::ParticleType t, double range, algebra::Vector3Ds sites) {
    sites_[t] = std::make_pair(range, sites);
  }
  static const char *get_name() {return "npctransport save";}
};

// for testing
IMPNPCTRANSPORTEXPORT void add_sites(RMF::FileHandle fh, core::ParticleType t,
                                     double radius, algebra::Vector3Ds sites);

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
                    atom::Hierarchy, atom::Hierarchies,
                    (RMF::FileConstHandle fh, Model *m),
                    See IMP::rmf::link_hierarchies() for more details.
                    The only difference is the addition of particle sites
                    support.
                    );

IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_RMF_LINKS_H */
