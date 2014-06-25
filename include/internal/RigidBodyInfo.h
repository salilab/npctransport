/**
 *  \file RigidBodyInfo.h
 *  \brief A summary of useful information about rigid bodies and their
 *         transformation for eg, caching purposes for SitesPairScore
 *
 *  Copyright 2007-2010 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_INTERNAL_RIGID_BODY_INFO_H
#define IMPNPCTRANSPORT_INTERNAL_RIGID_BODY_INFO_H

#include "../npctransport_config.h"
#include <IMP/core/rigid_bodies.h>
#include <IMP/base/check_macros.h>
#include <boost/current_function.hpp>

IMPNPCTRANSPORT_BEGIN_INTERNAL_NAMESPACE

#define INVALID_CACHE_ID 0

/** Reference frame information about rigid bodies and their
    transformation for eg, caching purposes for SitesPairScore.
    In particular saving the matrix representation of a roatation
    and its inverse are costly operations.
*/
struct RigidBodyInfo{
public:
  core::RigidBody rb;
  algebra::Transformation3D tr; // transformation of rb
  algebra::Rotation3D irot; //inverse rotation of rb
  unsigned int cache_id; // a serial number for keeping track of cache

public:
  //! updates the rigid body information from the stored rigid body
  //! and using passed cached id
  void update(unsigned int cache_id_arg)
  {
    cache_id = cache_id_arg;
    tr = rb.get_reference_frame().get_transformation_to();
    irot = tr.get_rotation().get_inverse();
  }

public:
  //! initiates an invalid rigid body info
RigidBodyInfo() : cache_id(INVALID_CACHE_ID) {}

  //! initiats the info for rigid body of particle index pi in model m,
  //! assigned with passed cached id. Can be overridden by set_particle
RigidBodyInfo(IMP::Model* m, ParticleIndex pi, unsigned int cache_id_arg)
: rb(m, pi)
  { update(cache_id_arg); }


  //! set/reset the info for rigid body of particle index pi in model m,
  //! assigned with passed cached id
  void set_particle(IMP::Model* m, ParticleIndex pi, unsigned int cache_id_arg)
  {
    IMP_USAGE_CHECK(core::RigidBody::get_is_setup(m, pi),
                    BOOST_CURRENT_FUNCTION << " pi must be a rigid body");
    rb = core::RigidBody(m,pi);
    update(cache_id_arg);
  }


};


IMPNPCTRANSPORT_END_INTERNAL_NAMESPACE

#endif /* IMPNPCTRANSPORT_INTERNAL_RIGID_BODY_INFO_H */
