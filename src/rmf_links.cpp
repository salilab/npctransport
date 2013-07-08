/**
 *  \file io.cpp
 *  \brief description.
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/rmf_links.h>
#include <IMP/core/rigid_bodies.h>
#include <IMP/npctransport/Transporting.h>
#include <IMP/core/XYZ.h>
#include <IMP/base/check_macros.h>
#include <IMP/base/exception.h>
#include <IMP/base/log.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE
HierarchyWithSitesLoadLink::HierarchyWithSitesLoadLink(RMF::FileConstHandle fh,
                                                       Model *m)
    : rmf::HierarchyLoadLink(fh, m), bf_(fh) {
  RMF::Category npc_cat = fh.get_category("npc");
  is_last_entry_from_top_key_ = fh.get_int_key(npc_cat, "last entry from top");
  n_entries_bottom_key_ = fh.get_int_key(npc_cat, "n entries bottom");
  n_entries_top_key_ = fh.get_int_key(npc_cat, "n entries top");
}

void HierarchyWithSitesLoadLink::do_load_node(RMF::NodeConstHandle nh,
                                              Particle *o) {
  rmf::HierarchyLoadLink::do_load_node(nh, o);
  // load particle transport directionality if needed
  if (nh.get_has_value(is_last_entry_from_top_key_)) {
    IMP_ALWAYS_CHECK(Transporting::get_is_setup(o),
                     "is_last_entry_from_top is relevant only for particles"
                     " decorated with Transporting class - particle " << *o,
                     IMP::base::ValueException);
    IMP_ALWAYS_CHECK(IMP::core::XYZ::get_is_setup(o) &&
                         nh.get_has_value(n_entries_bottom_key_) &&
                         nh.get_has_value(n_entries_top_key_),
                     "It is expected that a transporting particle would have "
                     "coordinates and the RMF node would have all relevant "
                     "keys for a transporting particle, particle " << *o,
                     IMP::base::ValueException);
    bool is_last_entry_from_top = nh.get_value(is_last_entry_from_top_key_);
    Transporting ot(o);
    ot.set_is_last_entry_from_top(is_last_entry_from_top);
    double cur_z = IMP::core::XYZ(o).get_coordinates()[2];
    ot.set_last_tracked_z(cur_z);
    ot.set_n_entries_bottom(nh.get_value(n_entries_bottom_key_));
    ot.set_n_entries_top(nh.get_value(n_entries_top_key_));
    IMP_LOG(VERBOSE, "Setting is_last_entry_from_top value to "
            << is_last_entry_from_top << "; last_tracked_z to " << cur_z
            << "; n_entries_bottom " << nh.get_value(n_entries_bottom_key_)
            << "; n_entries_top " << nh.get_value(n_entries_top_key_)
            << " - for particle " << *o << std::endl);
  }
}

void HierarchyWithSitesLoadLink::do_add_link_recursive(
    Particle *root, Particle *p, RMF::NodeConstHandle cur) {
  HierarchyLoadLink::do_add_link_recursive(root, p, cur);
  if (atom::Hierarchy(p).get_number_of_children() == 0 &&
      core::Typed::get_is_setup(p)) {
    core::ParticleType tp = core::Typed(p).get_type();
    algebra::Vector3Ds sites;
    RMF::NodeConstHandles children = cur.get_children();
    for (unsigned int i = 0; i < children.size(); ++i) {
      if (children[i].get_type() == RMF::GEOMETRY && bf_.get_is(children[i])) {
        RMF::BallConst b = bf_.get(children[i]);
        RMF::Floats cs = b.get_coordinates();
        sites.push_back(algebra::Vector3D(cs.begin(), cs.end()));
      }
    }
    if (sd_) sd_->set_sites(tp, sites);
  }
}

HierarchyWithSitesSaveLink::HierarchyWithSitesSaveLink(RMF::FileHandle fh)
    : rmf::HierarchySaveLink(fh), bf_(fh), cf_(fh) {
  RMF::Category npc_cat = fh.get_category("npc");
  is_last_entry_from_top_key_ = fh.get_int_key(npc_cat, "last entry from top");
  n_entries_bottom_key_ = fh.get_int_key(npc_cat, "n entries bottom");
  n_entries_top_key_ = fh.get_int_key(npc_cat, "n entries top");
}

std::pair<double, algebra::Vector3Ds> HierarchyWithSitesSaveLink::get_sites(
    core::ParticleType t) const {
  if (sd_) {
    return std::make_pair(sd_->get_site_radius(t), sd_->get_sites(t));
  } else if (sites_.find(t) != sites_.end()) {
    return sites_.find(t)->second;
  } else {
    return std::make_pair(0.0, algebra::Vector3Ds());
  }
}

void HierarchyWithSitesSaveLink::do_add_recursive(Particle *root, Particle *p,
                                                  RMF::NodeHandle cur) {
  if (!sd_ && p->has_attribute(get_simulation_data_key())) {
    Object *o = p->get_value(get_simulation_data_key());
    sd_ = dynamic_cast<SimulationData *>(o);
  }
  if (core::Typed::get_is_setup(p)) {
    core::ParticleType type = core::Typed(p).get_type();
    std::pair<double, algebra::Vector3Ds> sites = get_sites(type);
    for (unsigned int i = 0; i < sites.second.size(); ++i) {
      RMF::NodeHandle ch = cur.add_child("site", RMF::GEOMETRY);
      RMF::Floats color(3);
      color[0] = 1;
      color[1] = 0;
      color[2] = 0;
      cf_.get(ch).set_rgb_color(color);

      RMF::Ball b = bf_.get(ch);
      b.set_radius(sites.first);
      algebra::Vector3D local = sites.second[i];
      b.set_coordinates(
          RMF::Floats(local.coordinates_begin(), local.coordinates_end()));
    }
  }
  rmf::HierarchySaveLink::do_add_recursive(root, p, cur);
}

void HierarchyWithSitesSaveLink::do_save_node(Particle *p, RMF::NodeHandle n) {
  rmf::HierarchySaveLink::do_save_node(p, n);
  // update side of transport entry
  if (Transporting::get_is_setup(p)) {
    Transporting t(p);
    n.set_value(is_last_entry_from_top_key_, t.get_is_last_entry_from_top());
    n.set_value(n_entries_bottom_key_, t.get_n_entries_bottom());
    n.set_value(n_entries_top_key_, t.get_n_entries_top());
  }
}

IMP_DEFINE_LINKERS(HierarchyWithSites, hierarchy_with_sites,
                   hierarchies_with_sites, atom::Hierarchy, atom::Hierarchies,
                   atom::Hierarchy, atom::Hierarchies, (RMF::FileHandle fh),
                   (RMF::FileConstHandle fh, Model *m), (fh), (fh, m),
                   (fh, IMP::internal::get_model(hs)));

void add_sites(RMF::FileHandle fh, core::ParticleType t, double range,
               algebra::Vector3Ds sites) {
  HierarchyWithSitesSaveLink *l = get_hierarchy_with_sites_save_link(fh);
  l->add_sites(t, range, sites);
}

IMPNPCTRANSPORT_END_NAMESPACE
