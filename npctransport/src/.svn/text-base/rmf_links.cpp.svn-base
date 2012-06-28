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
HierarchyLoadLink::HierarchyLoadLink(RMF::FileConstHandle fh, Model *m):
  rmf::HierarchyLoadLink(fh, m){}

HierarchySaveLink::HierarchySaveLink(RMF::FileHandle fh):
  rmf::HierarchySaveLink(fh), bf_(fh), cf_(fh) {}

void HierarchySaveLink::do_add_recursive(Particle *root,
                                         Particle *p, RMF::NodeHandle cur) {
  if (core::Typed::particle_is_instance(p)) {
    core::ParticleType type= core::Typed(p).get_type();
    algebra::Vector3Ds sites= sd_->get_sites(type);
    for (unsigned int i=0; i< sites.size(); ++i) {
      cur.add_child("site", RMF::GEOMETRY);
    }
  }
  if (!sd_ && p->has_attribute(get_simulation_data_key())) {
    Object *o= p->get_value(get_simulation_data_key());
    sd_= dynamic_cast<SimulationData*>(o);
  }
  rmf::HierarchySaveLink::do_add_recursive(root, p, cur);
}
void HierarchySaveLink::do_save_node(Particle *p,
                                     RMF::NodeHandle n,
                                     unsigned int frame) {
  if (core::Typed::particle_is_instance(p)) {
    core::ParticleType type= core::Typed(p).get_type();
    algebra::Vector3Ds sites= sd_->get_sites(type);
    core::RigidBody rb(p);
    RMF::NodeHandles ch= n.get_children();
    algebra::ReferenceFrame3D rf= rb.get_reference_frame();
    for (unsigned int i=0; i< sites.size(); ++i) {
      RMF::Ball b= bf_.get(ch[i], frame);
      b.set_radius(sd_->get_site_radius(type));
      algebra::Vector3D local= rf.get_global_coordinates(sites[i]);
      b.set_coordinates(RMF::Floats(local.coordinates_begin(),
                                    local.coordinates_end()));
      RMF::Floats color(3);
      color[0]=1; color[1]=0; color[2]=0;
      cf_.get(ch[i], frame).set_rgb_color(color);
    }
  }
  rmf::HierarchySaveLink::do_save_node(p,n,frame);
}


IMP_DEFINE_LINKERS(Hierarchy, hierarchy, hierarchies,
                   atom::Hierarchy,atom::Hierarchies,
                   atom::Hierarchy,atom::Hierarchies,
                   (RMF::FileHandle fh),
                   (RMF::FileConstHandle fh,
                    Model *m), (fh), (fh, m),
                   (fh, IMP::internal::get_model(hs)));

IMPNPCTRANSPORT_END_NAMESPACE
