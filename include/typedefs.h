/**
 *  \file typedefs.h
 *  \brief description
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_TYPEDEFS_H
#define IMPNPCTRANSPORT_TYPEDEFS_H

#include "npctransport_config.h"
#include <IMP/core/Typed.h>
#include <utility>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

// an interaction that involves particles of two types
typedef std::pair
<IMP::core::ParticleType, IMP::core::ParticleType> InteractionType;


IMPNPCTRANSPORT_END_NAMESPACE

#endif  /* IMPNPCTRANSPORT_TYPEDEFS_H */
