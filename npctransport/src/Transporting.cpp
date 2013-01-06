/**
 *  \file Transporting.cpp
 *  \brief description.
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/Transporting.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

IntKey Transporting::get_is_last_entry_from_top_key() {
  static IntKey ik("last entry is from top");
  return ik;
}

FloatKey Transporting::get_last_tracked_z_key() {
  static FloatKey fk("last tracked z");
  return fk;
}

void Transporting::show(std::ostream &out) const
{
  out << "Transporting is_last_entry_from_top = "
      << get_is_last_entry_from_top() <<
    " ; last_tracked_z = "
      << get_last_tracked_z();
}

IMPNPCTRANSPORT_END_NAMESPACE
