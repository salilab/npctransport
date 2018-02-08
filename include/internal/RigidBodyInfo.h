/**
 *  \file RigidBodyInfo.h
 *  \brief A summary of useful information about rigid bodies and their
 *         transformation for eg, caching purposes for SitesPairScore
 *
 *  Copyright 2007-2018 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_INTERNAL_RIGID_BODY_INFO_H
#define IMPNPCTRANSPORT_INTERNAL_RIGID_BODY_INFO_H

#include "../npctransport_config.h"
#include <IMP/algebra/Sphere3D.h>
#include <IMP/core/rigid_bodies.h>
#include <IMP/check_macros.h>
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
  ParticleIndex pi;
  algebra::Transformation3D tr; // transformation of rb
  algebra::Rotation3D irot; //inverse rotation of rb
  double radius;
  double iradius; // inverse radius
  unsigned int cache_id; // a serial number for keeping track of cache

public:
  //! updates the rigid body information from the stored rigid body
  //! and using passed cached id
  //!
  //! @param st - the table of spheres data in model for each particle index,
  //!             e.g. st[pi.get_index()] has the sphere of particle pi in the
  //!             model of rb.
  //! @param qt - an array of arrays for each quaternion component,
  //!             e.g. qt[i][pi.get_index()] has the i'th quaternion component
  //!             of particle pi in the model of rb
  //! @param cache_id_arg - the new cache id asociated with this rigid body
  void update
  (algebra::Sphere3D const* st, // spheres table
   double const **qt, // quaternions table
   unsigned int cache_id_arg);

public:
  //! initiates an invalid rigid body info
RigidBodyInfo() : cache_id(INVALID_CACHE_ID) {}

  //! initiats the info for rigid body of particle index pi in model m,
  //! assigned with passed cached id. Can be overridden by set_particle
RigidBodyInfo(algebra::Sphere3D const* spheres_table,
              double const**quaterntions_tables,
              ParticleIndex pi,
              unsigned int cache_id_arg)
: pi(pi)
  {
    update(spheres_table, quaterntions_tables, cache_id_arg);
  }

  //! set/reset the info for rigid body of particle index pi in model m,
  //! assigned with passed cached id
  void set_particle
  (algebra::Sphere3D const* spheres_table,
   double const**quaternions_tables,
   ParticleIndex pi,
   unsigned int cache_id_arg)
  {
    //    IMP_USAGE_CHECK(core::RigidBody::get_is_setup(m, pi),
    //                BOOST_CURRENT_FUNCTION << " pi must be a rigid body");
    this->pi=pi;
    update(spheres_table, quaternions_tables, cache_id_arg);
  }


};

//!
inline
void RigidBodyInfo::update
(algebra::Sphere3D const* st, // spheres table
   double const **qt, // quaternions table
   unsigned int cache_id_arg)
  {
    this->cache_id = cache_id_arg;
    algebra::Sphere3D s=st[pi.get_index()];
    algebra::Rotation3D rot(qt[0][pi.get_index()],
                            qt[1][pi.get_index()],
                            qt[2][pi.get_index()],
                            qt[3][pi.get_index()]);
    this->tr=algebra::Transformation3D(rot, s.get_center());
    this->irot = rot.get_inverse();
    this->radius=s.get_radius();
    this->iradius = 1.0/radius;
  }


IMPNPCTRANSPORT_END_INTERNAL_NAMESPACE

#endif /* IMPNPCTRANSPORT_INTERNAL_RIGID_BODY_INFO_H */
