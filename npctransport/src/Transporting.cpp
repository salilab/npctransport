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


void Transporting::show(std::ostream &out) const
{
  out << "Transporting is_last_entry_from_top = "
      << get_is_last_entry_from_top();
}

IMPNPCTRANSPORT_END_NAMESPACE
