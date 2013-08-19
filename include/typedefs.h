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
#include <iostream>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

// an interaction that involves particles of two types
typedef std::pair<IMP::core::ParticleType, IMP::core::ParticleType>
    InteractionType;

typedef IMP::base::set< core::ParticleType > ParticleTypeSet ;

IMPNPCTRANSPORT_END_NAMESPACE

#ifndef SWIG
inline std::ostream& operator<<
( std::ostream& o, IMP::npctransport::InteractionType const& it );

inline std::ostream& operator<<
( std::ostream& o, IMP::npctransport::InteractionType const& it )
{
  o << "(" << it.first << "," << it.second << ")";
  return o;
}
#endif


#endif /* IMPNPCTRANSPORT_TYPEDEFS_H */
