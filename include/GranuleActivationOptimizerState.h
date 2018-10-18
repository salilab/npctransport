/**
 * \file pbc/GranuleActivationOptimizerState.h
 * \brief an optimizer state that rigidifies the glucose and insulin granules upon collision, also,
 *  deactivate the bound binding sites on insulin granules.
 *
 * Description:
 * 1. Get optimizer state for each frame of the trajectory (glucose and granule diffusion).
 * 2. Glucose form a rigid body with a granule[i] and its patches, upon collision with patches.
 * 3. The bound patches are deactivated.
 * 4. Update the optimizer state.
 *
 *
 *  Copyright 2007-2018 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_GRANULE_ACTIVATION_OPTIMIZER_STATE_H
#define IMPNPCTRANSPORT_GRANULE_ACTIVATION_OPTIMIZER_STATE_H

#include "npctransport_config.h"
#include <IMP/OptimizerState.h>
#include <IMP/SingletonContainer.h>
#include <IMP/container/CloseBipartitePairContainer.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE
//TODO: add tester

/**
 */
class IMPNPCTRANSPORTEXPORT GranuleActivationOptimizerState
    : public OptimizerState {
 private:
  typedef OptimizerState P;

  // maintains a list of nearby particle pairs in a bipartite graph
  IMP::PointerMember<IMP::container::CloseBipartitePairContainer>
      close_bipartite_pair_container_;

  int periodicity_;

 public:
  /**
     An optimizer states that rigidifies nearby vesicles and glucose molecules if certain conditions are met (e.g. interaction)

     @param vesicles_container container of diffusing vesicles (which may change dynamically after construction)
     @param glucose_container container of diffusing glucose molecules (which may change dynamically after construction)
     @param contact_range the range of sphere distance in angstroms under which vesicles and glucose will be tested for rigidification
     @param slack slack in angstroms for an IMP::container::CloseBipartitePairContainer object that tracks nearby vesicles and glucose
                  (affects running time - large slack may result in a longer list of pairs being tracked, up to contact_range+2*slack,
                  but the maintenance of CloseBipartitePairContainer may be faster)
     @param periodicity the frame interval for updating this optimizer state
   */
  GranuleActivationOptimizerState
    ( IMP::SingletonContainerAdaptor vesicles_container,
      IMP::SingletonContainerAdaptor glucose_container,
      double contact_range,
      double slack = 1.0,
      unsigned int periodicity=1);

 protected:
  virtual void do_update(unsigned int call_num) IMP_OVERRIDE;

 private:
  /**
     get the index of the interacting patch of pip[0] with pip[1],
     if they are interacting, or -1 if they are not interacting

     @param pip an ordered pair of particle indexes in the model
            of the optimizer state, pip[0] is assumed to be
            a vesicle particle and pip[1] is a glucose molecule
   */
  int get_interacting_patch(ParticleIndexPair pip) const;

  /*
     @param pip an ordered pair of particle indexes in the model
            of the optimizer state, pip[0] is assumed to be
            a vesicle particle and pip[1] is a glucose molecule.
  */
  void rigidify_pair(ParticleIndexPair pip);

 public:
  IMP_OBJECT_METHODS(GranuleActivationOptimizerState);
};
IMP_OBJECTS(GranuleActivationOptimizerState, GranuleActivationOptimizerStates);

IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_GRANULE_ACTIVATION_OPTIMIZER_STATE_H */
