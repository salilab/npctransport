/**
 *  \file GranuleActivationOptimizerState.cpp
 *  \brief
 *
 *  Copyright 2007-2018 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/GranuleActivationOptimizerState.h>
#include <limits>

IMPNPCTRANSPORT_BEGIN_NAMESPACE
GranuleActivationOptimizerState::GranuleActivationOptimizerState
( IMP::SingletonContainer* vesicles_container,
  IMP::SingletonContainer::SingletonContainer* glucose_container,
  double contact_range, double slack = 1.0,
  unsigned int periodicity=1);
: P(vesicles_container ? vesicles_container->get_model() :  nullptr,
      "GranuleActivationOptimizerState%1%"),
{
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
  IMP_UNUSED(call_num);
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
  // TODO: right now fake interaction with first patch, need to implement actual detection
  //       of interacting patch
  return 0;
}


void
GranuleActivationOptimizerState::rigidify_pair
(ParticleIndexPair pip)
{
  Model* m= get_model();
  if(!core::RigidBody::get_is_setup(m, pip[0])){
    core::RigidBody::do_setup_particle(m, pip[0]);
  }
  core::RigidBody rb0= core::RigidBody(m, pip[0]);
  rb0.add_member(pip[1])
}

IMPNPCTRANSPORT_END_NAMESPACE
