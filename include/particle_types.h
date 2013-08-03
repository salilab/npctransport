/**
 *  \file particle_types.h
 *  \brief description
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_PARTICLE_TYPES_H
#define IMPNPCTRANSPORT_PARTICLE_TYPES_H

#include "npctransport_config.h"
#include <IMP/atom/Hierarchy.h>
#include <IMP/core/Typed.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

IMPNPCTRANSPORT_DEPRECATED_HEADER
(2.1, "header is only left for backward comaptability with old inputs - types are now specified in protobuf")

IMPNPCTRANSPORTEXPORT
extern core::ParticleTypes type_of_fg;
IMPNPCTRANSPORTEXPORT
extern core::ParticleTypes type_of_float;

inline
IMPNPCTRANSPORT_DEPRECATED_FUNCTION_DECL(2.1)
 IMPNPCTRANSPORTEXPORT
core::ParticleType get_type_of_fg(
    unsigned int index) {
  IMPNPCTRANSPORT_DEPRECATED_FUNCTION_DEF(2.1, "all in protobuf now");
  return type_of_fg[index];
}

inline
IMPNPCTRANSPORT_DEPRECATED_FUNCTION_DECL(2.1)
 IMPNPCTRANSPORTEXPORT
unsigned int get_number_of_types_of_fg() {
  IMPNPCTRANSPORT_DEPRECATED_FUNCTION_DEF(2.1, "all in protobuf now");
  return type_of_fg.size();
}

inline
IMPNPCTRANSPORT_DEPRECATED_FUNCTION_DECL(2.1)
IMPNPCTRANSPORTEXPORT core::ParticleType get_type_of_float(
    unsigned int index) {
  IMPNPCTRANSPORT_DEPRECATED_FUNCTION_DEF(2.1, "all in protobuf now");
  return type_of_float[index];
}

inline
IMPNPCTRANSPORT_DEPRECATED_FUNCTION_DECL(2.1)
IMPNPCTRANSPORTEXPORT unsigned int get_number_of_types_of_float() {
  IMPNPCTRANSPORT_DEPRECATED_FUNCTION_DEF(2.1, "all in protobuf now");
  return type_of_float.size();
}

/** \deprecated_at{2.1} Use the SimulationData::get_fg_chains() method
    instead
*/
IMPNPCTRANSPORT_DEPRECATED_FUNCTION_DECL(2.1)
IMPNPCTRANSPORTEXPORT
atom::Hierarchies get_fg_chains(atom::Hierarchy h);

IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_PARTICLE_TYPES_H */
