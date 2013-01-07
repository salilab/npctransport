/**
 *  \file Transporting.cpp
 *  \brief description.
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/Transporting.h>
#include <IMP/core/XYZ.h>
#include <IMP/base/exception.h>
#include <IMP/base/check_macros.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

Transporting Transporting::setup_particle(Particle *p,
                                          bool is_last_entry_from_top)
{
  int hinum = 0;
  IMP_ALWAYS_CHECK(IMP::core::XYZ::particle_is_instance(p),
                   "It is expected that a transporting particle would have "
                   "coordinates, particle " << *p,
                   IMP::base::ValueException);
  p->add_attribute(get_is_last_entry_from_top_key(), is_last_entry_from_top);
  double cur_z = IMP::core::XYZ(p).get_coordinates()[2];
  p->add_attribute(get_last_tracked_z_key(), cur_z);
  p->add_attribute(get_n_entries_bottom_key(), 0);
  p->add_attribute(get_n_entries_top_key(), 0);
  return Transporting(p);
}

IntKey Transporting::get_is_last_entry_from_top_key() {
  static IntKey ik("last entry is from top");
  return ik;
}

FloatKey Transporting::get_last_tracked_z_key() {
  static FloatKey fk("last tracked z");
  return fk;
}

IntKey Transporting::get_n_entries_bottom_key() {
  static IntKey ik("n_entries_bottom");
  return ik;
}

IntKey Transporting::get_n_entries_top_key() {
  static IntKey ik("n_entries_top");
  return ik;
}

void Transporting::show(std::ostream &out) const
{
  out << "Transporting is_last_entry_from_top = "
      << get_is_last_entry_from_top() <<
    " ; last_tracked_z = "
      << get_last_tracked_z();
}

IMPNPCTRANSPORT_END_NAMESPACE
