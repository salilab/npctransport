/**
 *  \file Transporting.cpp
 *  \brief description.
 *
 *  Copyright 2007-2017 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/Transporting.h>
#include <IMP/core/XYZ.h>
#include <IMP/exception.h>
#include <IMP/check_macros.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

void Transporting::do_setup_particle(IMP::Model* m,
                                             ParticleIndex pi,
                                             bool is_last_entry_from_top) {
  IMP_ALWAYS_CHECK(IMP::core::XYZ::get_is_setup(m, pi),
                   "It is expected that a transporting particle would have "
                   "coordinates, particle index " << pi,
                   IMP::ValueException);
  double cur_z = IMP::core::XYZ(m,pi).get_coordinates()[2];
  m->add_attribute(get_last_tracked_z_key(), pi, cur_z);
  m->add_attribute(get_n_entries_bottom_key(), pi, 0);
  m->add_attribute(get_n_entries_top_key(), pi, 0);
  m->add_attribute(get_is_last_entry_from_top_key(), pi, is_last_entry_from_top);
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

void Transporting::show(std::ostream &out) const {
  out << "Transporting is_last_entry_from_top = "
      << get_is_last_entry_from_top()
      << " ; last_tracked_z = " << get_last_tracked_z();
}

IMPNPCTRANSPORT_END_NAMESPACE
