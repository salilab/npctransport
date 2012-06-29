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


IMPNPCTRANSPORTEXPORT
extern core::ParticleTypes type_of_fg;
IMPNPCTRANSPORTEXPORT
extern core::ParticleTypes type_of_float;

inline IMPNPCTRANSPORTEXPORT
core::ParticleType get_type_of_fg(unsigned int index)
{
  return type_of_fg[index];
}

inline IMPNPCTRANSPORTEXPORT
core::ParticleType get_type_of_float(unsigned int index)
{
  return type_of_float[index];
}

IMPNPCTRANSPORTEXPORT
atom::Hierarchies get_fg_chains(atom::Hierarchy h);

IMPNPCTRANSPORT_END_NAMESPACE

#endif  /* IMPNPCTRANSPORT_PARTICLE_TYPES_H */
