/**
 *  \file GranuleActivationOptimizerState.cpp
 *  \brief
 *
 *  Copyright 2007-2018 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/GranuleActivationOptimizerState.h>
#include <IMP/algebra/Transformation3D.h>
#include <IMP/algebra/ReferenceFrame3D.h>
#include <IMP/core/rigid_bodies.h>
#include <limits>

IMPNPCTRANSPORT_BEGIN_NAMESPACE
GranuleActivationOptimizerState::GranuleActivationOptimizerState
( IMP::SingletonContainer* vesicles_container,
  IMP::SingletonContainer* glucose_container,
  double contact_range,
  double slack,
  unsigned int periodicity)
: P(vesicles_container ? vesicles_container->get_model() :  nullptr,
    "GranuleActivationOptimizerState%1%"),
  periodicity_(periodicity)
{
  IMP_UNUSED(periodicity);
  close_bipartite_pair_container_ =
    new IMP::container::CloseBipartitePairContainer
    (vesicles_container,
     glucose_container,
     contact_range,
     slack);
}

void GranuleActivationOptimizerState::do_update
(unsigned int call_num)
{
  if(call_num % periodicity_ != 0){
    return;
  }
  IMP_CONTAINER_FOREACH                                 \
    (IMP::container::CloseBipartitePairContainer,
     close_bipartite_pair_container_,
     {
       ParticleIndexPair const& pip = _1;
       int ipatch= get_interacting_patch(pip);
       if(ipatch>=0){
         rigidify_pair(pip);
         // TODO: need to also deactivate patch ipatch
       }
     }
     );
}

int
GranuleActivationOptimizerState::get_interacting_patch
(ParticleIndexPair pip) const
{
  IMP_UNUSED(pip);
  // TODO: right now fake interaction with first patch, need to implement actual detection
  //       of interacting patch
  return 0;
}


void
GranuleActivationOptimizerState::rigidify_pair
(ParticleIndexPair pip)
{
  Model* m= get_model();
  IMP_USAGE_CHECK(core::XYZR::get_is_setup(m, pip[0]),
                  "particles for rigidifications must be spheres as well");
  core::XYZR xyzr0(m, pip[0]);
  if(!core::RigidBody::get_is_setup(m, pip[0])){
    algebra::Transformation3D tr(xyzr0.get_sphere().get_center());
    algebra::ReferenceFrame3D rf(tr);
    core::RigidBody::setup_particle(m, pip[0], rf);
  }
  // Add rigid body member but make sure to keep original radius for this particle:
  core::RigidBody rb0= core::RigidBody(m, pip[0]);
  double original_radius= xyzr0.get_radius();
  rb0.add_member(pip[1]);
  xyzr0.set_radius(original_radius);
  xyzr0.set_coordinates_are_optimized(false);
}

IMPNPCTRANSPORT_END_NAMESPACE
