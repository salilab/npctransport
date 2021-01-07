/**
 *  \file typedefs.h
 *  \brief description
 *
 *  Copyright 2007-2021 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_TYPEDEFS_H
#define IMPNPCTRANSPORT_TYPEDEFS_H

#include "npctransport_config.h"
#include <IMP/core/Typed.h>
#include <boost/unordered_set.hpp>
#include <utility>
#include <iostream>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

//! an interaction that involves particles of two types
typedef std::pair<IMP::core::ParticleType, IMP::core::ParticleType>
    InteractionType;
IMP_VALUES(InteractionType, InteractionTypes);

//! convenience method for creating and interaction type in swig
inline IMPNPCTRANSPORTEXPORT
InteractionType make_ordered_interaction_type(IMP::core::ParticleType t0,
                                              IMP::core::ParticleType t1)
{
  return InteractionType(t0,t1);
}

//! an interaction type canonized to (t0,t1) s.t. t0 >= t1,
//! so order doesn't matter
inline IMPNPCTRANSPORTEXPORT
InteractionType make_unordered_interaction_type(IMP::core::ParticleType t0,
                                                IMP::core::ParticleType t1)
{
  return (t0 >= t1) ? InteractionType(t0,t1) : InteractionType(t1,t0);
}


//! ParticleIndexPair canonized to (pi0',pi1') s.t. pi0' >= pi1'
//! so order doesn't matter
inline IMP::ParticleIndexPair make_unordered_particle_index_pair
(IMP::ParticleIndex pi0, IMP::ParticleIndex pi1)
{
  return (pi0 >= pi1)
    ? ParticleIndexPair(pi0, pi1)
    : ParticleIndexPair(pi1, pi0);
}

//! ParticleIndexPair canonized to (pi0',pi1') s.t. pi0' >= pi1'
//! so order doesn't matter
inline ParticleIndexPair make_unordered_particle_index_pair
(ParticleIndexPair pip)
{
  return make_unordered_particle_index_pair(pip[0], pip[1]);
}

typedef boost::unordered_set<core::ParticleType> ParticleTypeSet;

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
