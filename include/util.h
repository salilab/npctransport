/**
# * \file util.h
# * \brief utility methods for npctransport

Simulate an fg and a kap interacting
#
# * Copyright 2007-2012 IMP Inventors. All rights reserved.
# */

#ifndef IMPNPCTRANSPORT_UTIL_H
#define IMPNPCTRANSPORT_UTIL_H

#include "npctransport_config.h"
#include <IMP/npctransport/npctransport_proto.fwd.h>
#include <IMP/core/Typed.h>
#include <IMP/core/rigid_bodies.h>
#include <IMP/algebra/Vector3D.h>
#include <IMP/algebra/ReferenceFrame3D.h>
#include <IMP/algebra/Transformation3D.h>
#include "typedefs.h"
#include <IMP/base_types.h>
#include <string>

IMPNPCTRANSPORT_BEGIN_NAMESPACE


/** returns particles with optimizable coordinates from particles */
ParticlesTemp get_optimizable_particles
(ParticlesTemp const& particles);


#ifndef SWIG

/** finds the index of s.fgs() whose type equals pt.get_string()
    if it does not exist, add it to s
    @param s the statistics message for searching the fg
    @param pt the type of fg to look for
*/
unsigned int find_or_add_fg_of_type(::npctransport_proto::Statistics* s,
                                    IMP::core::ParticleType pt);

/** finds the index of s->floaters() whose type equals pt.get_string()
    if it does not exist, add it to s
    @param s the statistics message for searching t
    @param pt the type to look for
*/
unsigned int find_or_add_floater_of_type(::npctransport_proto::Statistics* s,
                                         IMP::core::ParticleType pt);


/** finds the index of s.interactions() whose type0 and type1 particle
  type equivalents are equale to it. If it does not exist, add it to s
  @param s the statistics message for searching t
  @param it the type to look for
  @return the index of the interaction type in s.interaction()
*/
unsigned int find_or_add_interaction_of_type
( ::npctransport_proto::Statistics* s,
  IMP::npctransport::InteractionType it);

/**
   @param p a rigid body particle
   @param local the vector in local coordinates

   @return the global coordinates of local based on the reference frame of the
           rigid body p
*/
inline algebra::Vector3D get_global_from_local_v3( Particle* p,
                                            const algebra::Vector3D& local);


algebra::Vector3D get_global_from_local_v3( Particle* p,
                                              const algebra::Vector3D& local)
{
  IMP_USAGE_CHECK(core::RigidBody::particle_is_instance(p),
                  "Particle must be rigid body in order to use"
                    " its ref frame in get_global_from_local_v3");
  core::RigidBody rb(p);
  algebra::ReferenceFrame3D rf = rb.get_reference_frame();
  algebra::Transformation3D t = rf.get_transformation_to();
  return t.get_transformed(local);
}

#endif  // SWIG

IMPNPCTRANSPORT_END_NAMESPACE


#endif /* IMPNPCTRANSPORT_UTIL_H */
