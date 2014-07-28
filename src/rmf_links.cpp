/**
 *  \file io.cpp
 *  \brief description.
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/rmf_links.h>
#include <IMP/npctransport/util.h>

#include <IMP/atom/Hierarchy.h>
#include <IMP/core/rigid_bodies.h>
#include <IMP/npctransport/Transporting.h>
#include <IMP/core/XYZ.h>
#include <IMP/base/check_macros.h>
#include <IMP/base/exception.h>
#include <IMP/base/log.h>

#include <boost/foreach.hpp>

IMPNPCTRANSPORT_BEGIN_NAMESPACE
HierarchyWithSitesLoadLink::HierarchyWithSitesLoadLink(RMF::FileConstHandle fh)
    : rmf::HierarchyLoadLink(fh), bf_(fh) {
  RMF::Category npc_cat = fh.get_category("npc");
  is_last_entry_from_top_key_ =
      fh.get_key<RMF::IntTraits>(npc_cat, "last entry from top");
  n_entries_bottom_key_ =
      fh.get_key<RMF::IntTraits>(npc_cat, "n entries bottom");
  n_entries_top_key_ = fh.get_key<RMF::IntTraits>(npc_cat, "n entries top");
  RMF::Category imp_cat = fh.get_category("imp");
  coordinates_are_optimized_key_ =
    fh.get_key<RMF::IntTraits>(imp_cat,"coordinates_are_optimized");
}

void HierarchyWithSitesLoadLink::do_load_hierarchy(
    RMF::NodeConstHandle root_node, Model *m,
    ParticleIndex root) {
  BOOST_FOREACH(ParticleIndex pi, particles_.find(root)->second) {
    // load particle transport directionality if needed
    RMF::NodeConstHandle nh = rmf::get_node_from_association(
        root_node.get_file(), m->get_particle(pi));
    // TRANSPORTING attributes:
    if (nh.get_has_value(is_last_entry_from_top_key_)) {
      IMP_ALWAYS_CHECK(Transporting::get_is_setup(m, pi),
                       "is_last_entry_from_top is relevant only for particles"
                       " decorated with Transporting class - particle "
                       << m->get_particle_name(pi),
                       IMP::base::ValueException);
      IMP_ALWAYS_CHECK(IMP::core::XYZ::get_is_setup(m, pi) &&
                       nh.get_has_value(n_entries_bottom_key_) &&
                       nh.get_has_value(n_entries_top_key_),
                       "It is expected that a transporting particle would have "
                       "coordinates and the RMF node would have all relevant "
                       "keys for a transporting particle, particle "
                       << m->get_particle_name(pi),
                       IMP::base::ValueException);
      bool is_last_entry_from_top = nh.get_value(is_last_entry_from_top_key_);
      Transporting ot(m, pi);
      ot.set_is_last_entry_from_top(is_last_entry_from_top);
      double cur_z = IMP::core::XYZ(m, pi).get_coordinates()[2];
      ot.set_last_tracked_z(cur_z);
      ot.set_n_entries_bottom(nh.get_value(n_entries_bottom_key_));
      ot.set_n_entries_top(nh.get_value(n_entries_top_key_));
      IMP_LOG(VERBOSE, "Setting is_last_entry_from_top value to "
              << is_last_entry_from_top << "; last_tracked_z to " << cur_z
              << "; n_entries_bottom " << nh.get_value(n_entries_bottom_key_)
              << "; n_entries_top " << nh.get_value(n_entries_top_key_)
              << " - for particle " << m->get_particle_name(pi) << std::endl);
    }
    // IS OPTIMIZED attribute:
    if (nh.get_has_value(coordinates_are_optimized_key_)) {
      IMP_ALWAYS_CHECK(core::XYZ::get_is_setup(m, pi),
                       "coordinates are optimized are only valid"
                       " for XYZ decorated particles",
                       IMP::base::ValueException);
      core::XYZ xyz(m, pi);
      xyz.set_coordinates_are_optimized
        ( nh.get_static_value( coordinates_are_optimized_key_ ) );
    }
  }
}

void HierarchyWithSitesLoadLink::do_link_particle(Model *m,
                                                  ParticleIndex root,
                                                  ParticleIndex p,
                                                  RMF::NodeConstHandle node) {
  if (atom::Hierarchy(m, p).get_number_of_children() == 0 &&
      core::Typed::get_is_setup(m, p)) {
    particles_[root].push_back(p);
    core::ParticleType tp = core::Typed(m, p).get_type();
    algebra::Sphere3Ds sites;
    RMF::NodeConstHandles children = node.get_children();
    for (unsigned int i = 0; i < children.size(); ++i) {
      if (children[i].get_type() == RMF::GEOMETRY && bf_.get_is(children[i])) {
        RMF::decorator::BallConst b = bf_.get(children[i]);
        RMF::Vector3 cs = b.get_coordinates();
        algebra::Vector3D v(cs.begin(), cs.end());
        sites.push_back(algebra::Sphere3D(v, 0.0)); // TODO: use stored radius?
      }
    }
    if (sd_ && sites.size() > 0){
      sd_->set_sites(tp, sites);
    }
  }
}

HierarchyWithSitesSaveLink::HierarchyWithSitesSaveLink(RMF::FileHandle fh)
    : rmf::HierarchySaveLink(fh), bf_(fh), cf_(fh) {
  RMF::Category npc_cat = fh.get_category("npc");
  is_last_entry_from_top_key_ =
    fh.get_key<RMF::IntTraits>(npc_cat, "last entry from top");
  n_entries_bottom_key_ =
      fh.get_key<RMF::IntTraits>(npc_cat, "n entries bottom");
  n_entries_top_key_ = fh.get_key<RMF::IntTraits>(npc_cat, "n entries top");
  RMF::Category imp_cat = fh.get_category("imp");
  coordinates_are_optimized_key_ =
    fh.get_key<RMF::IntTraits>(imp_cat,"coordinates_are_optimized");
}

void HierarchyWithSitesSaveLink::add_sites_to_node
( RMF::NodeHandle cur_node, core::ParticleType t) const
{
  algebra::Sphere3Ds S;
  if (sd_)
    {
      S = sd_->get_sites(t);
    }
  else if (test_sites_.find(t) != test_sites_.end())
    { // note - sites_ is only for external testing w/o sd
      S = test_sites_.find(t)->second;
    }
  for (unsigned int i = 0; i < S.size(); ++i) {
    RMF::NodeHandle ch = cur_node.add_child("site", RMF::GEOMETRY);
    RMF::decorator::Ball b = bf_.get(ch);
    b.set_coordinates(RMF::Vector3(S[i].get_center()));
    double r = S[i].get_radius();
    if(r <= 0.0){
      // TODO: always use real radius? problematic cause 0.0 many times - just cause of display in chimera
      //       maybe find a solution with chimera - or a different field for backward compatibility
      if(sd_){
        r = sd_->get_site_display_radius(t);
      } else {
        r = 1.0;
      }
    }
    b.set_radius(r);
    RMF::Vector3 color(1,0,0);
    cf_.get(ch).set_rgb_color(color);
  }
}

void HierarchyWithSitesSaveLink::do_setup_node(kernel::Model *m,
                                             kernel::ParticleIndex root,
                                             kernel::ParticleIndex cur,
                                             RMF::NodeHandle cur_node) {
  if (!sd_ && m->get_has_attribute(get_simulation_data_key(), cur)) {
    Object *o = m->get_attribute(get_simulation_data_key(), cur);
    sd_ = dynamic_cast<SimulationData *>(o);
  }
  if (core::Typed::get_is_setup(m, cur)) {
    particles_[root].push_back(cur);
    IMP_USAGE_CHECK(atom::Hierarchy::get_is_setup(m, cur),
                    cur << " is expected to by Hierarchy decorated");
    if( atom::Hierarchy(m, cur).get_number_of_children() == 0 ) // leaf
      { // TODO: test leafiness / beadiness through sim-data?
        core::ParticleType type = core::Typed(m, cur).get_type();
        add_sites_to_node(cur_node, type);
      }
  }
}

void HierarchyWithSitesSaveLink::do_save_hierarchy(Model *m,
                                                   ParticleIndex root,
                                                   RMF::NodeHandle root_node) {
  BOOST_FOREACH(ParticleIndex p, particles_.find(root)->second) {
    RMF::NodeHandle n = rmf::get_node_from_association(root_node.get_file(),
                                                 m->get_particle(p));
    if (Transporting::get_is_setup(m, p)) {
      Transporting t(m, p);
      n.set_value(is_last_entry_from_top_key_, t.get_is_last_entry_from_top());
      n.set_value(n_entries_bottom_key_, t.get_n_entries_bottom());
      n.set_value(n_entries_top_key_, t.get_n_entries_top());
    }
    if (core::XYZ::get_is_setup(m, p)) {
      core::XYZ xyz(m, p);
      n.set_static_value(coordinates_are_optimized_key_,
                         xyz.get_coordinates_are_optimized());
    }
  }
}

IMP_DEFINE_LINKERS(HierarchyWithSites, hierarchy_with_sites,
                   hierarchies_with_sites, atom::Hierarchy, atom::Hierarchies,
                   (RMF::FileConstHandle fh, Model *m), (fh, m));

void add_test_sites(RMF::FileHandle fh, core::ParticleType t, double range,
               algebra::Vector3Ds sites) {
  HierarchyWithSitesSaveLink *l =
    rmf::internal::get_save_link<HierarchyWithSitesSaveLink>(fh);
  if(sites.size() > 0){
    l->add_test_sites(t, vectors2spheres(sites, range));
  }
}

//! for testing - adds the list of sites with specified display radius, to be
//! associated with particle type t. The file handle fh relies on this list
//! only if it doesn't have particles with simulation data keys
IMPNPCTRANSPORTEXPORT void add_test_sites(RMF::FileHandle fh,
                                          core::ParticleType t,
                                          algebra::Sphere3Ds sites)
 {
  HierarchyWithSitesSaveLink *l =
      rmf::internal::get_save_link<HierarchyWithSitesSaveLink>(fh);
  if(sites.size() > 0){
    l->add_test_sites(t, sites);
  }
}


IMPNPCTRANSPORT_END_NAMESPACE
