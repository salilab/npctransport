/**
 *  \file io.cpp
 *  \brief description.
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/rmf_links.h>
#include <IMP/core/rigid_bodies.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE
HierarchyWithSitesLoadLink::HierarchyWithSitesLoadLink(RMF::FileConstHandle fh, Model *m):
  rmf::HierarchyLoadLink(fh, m){}

HierarchyWithSitesSaveLink::HierarchyWithSitesSaveLink(RMF::FileHandle fh):
  rmf::HierarchySaveLink(fh), bf_(fh), cf_(fh) {}

std::pair<double, algebra::Vector3Ds>
HierarchyWithSitesSaveLink::get_sites(core::ParticleType t) const  {
  if (sd_) {
    return std::make_pair(sd_->get_site_radius(t),
                          sd_->get_sites(t) );
  } else if (sites_.find(t) != sites_.end()) {
    return sites_.find(t)->second;
  } else {
    return std::make_pair(0.0, algebra::Vector3Ds());
  }
}

void HierarchyWithSitesSaveLink::do_add_recursive(Particle *root,
                                         Particle *p, RMF::NodeHandle cur) {
  if (core::Typed::particle_is_instance(p)) {
    core::ParticleType type= core::Typed(p).get_type();
    unsigned int nsites= get_sites(type).second.size();
    for (unsigned int i=0; i< nsites; ++i) {
      cur.add_child("site", RMF::GEOMETRY);
    }
  }
  if (!sd_ && p->has_attribute(get_simulation_data_key())) {
    Object *o= p->get_value(get_simulation_data_key());
    sd_= dynamic_cast<SimulationData*>(o);
  }
  rmf::HierarchySaveLink::do_add_recursive(root, p, cur);
}
void HierarchyWithSitesSaveLink::do_save_node(Particle *p,
                                     RMF::NodeHandle n,
                                     unsigned int frame) {
  if (core::Typed::particle_is_instance(p)) {
    core::ParticleType type= core::Typed(p).get_type();
    std::pair<double, algebra::Vector3Ds> sites
        = get_sites(type);
    core::RigidBody rb(p);
    RMF::NodeHandles ch= n.get_children();
    algebra::ReferenceFrame3D rf= rb.get_reference_frame();
    for (unsigned int i=0; i< sites.second.size(); ++i) {
      RMF::Ball b= bf_.get(ch[i], frame);
      b.set_radius(sites.first);
      algebra::Vector3D local= rf.get_global_coordinates(sites.second[i]);
      b.set_coordinates(RMF::Floats(local.coordinates_begin(),
                                    local.coordinates_end()));
      RMF::Floats color(3);
      color[0]=1; color[1]=0; color[2]=0;
      cf_.get(ch[i], frame).set_rgb_color(color);
    }
  }
  rmf::HierarchySaveLink::do_save_node(p,n,frame);
}



IMP_DEFINE_LINKERS(HierarchyWithSites,
                   hierarchy_with_sites, hierarchies_with_sites,
                   atom::Hierarchy,atom::Hierarchies,
                   atom::Hierarchy,atom::Hierarchies,
                   (RMF::FileHandle fh),
                   (RMF::FileConstHandle fh, Model *m),
                   (fh), (fh, m),
                   (fh, IMP::internal::get_model(hs)));

void add_sites(RMF::FileHandle fh,
               core::ParticleType t,
               double range,
               algebra::Vector3Ds sites) {
  HierarchyWithSitesSaveLink *l=get_hierarchy_with_sites_save_link(fh);
  l->add_sites(t, range, sites);
}


IMPNPCTRANSPORT_END_NAMESPACE
