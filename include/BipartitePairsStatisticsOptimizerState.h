/**
 *  \file npctransport/BipartitePairsStatisticsOptimizerState.h
 *  \brief description
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_BIPARTITE_PAIRS_STATISTICS_OPTIMIZER_STATE_H
#define IMPNPCTRANSPORT_BIPARTITE_PAIRS_STATISTICS_OPTIMIZER_STATE_H

#include "npctransport_config.h"
#include <IMP/Particle.h>
#include <IMP/OptimizerState.h>
#include <IMP/optimizer_state_macros.h>
#include <IMP/core/PeriodicOptimizerState.h>
#include <IMP/container/CloseBipartitePairContainer.h>
#include <IMP/npctransport/typedefs.h>
#include <deque>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

/** Track the interaction between pairs from one group of particles
    with particles from another group, within some specified contat
    range
*/
class IMPNPCTRANSPORTEXPORT BipartitePairsStatisticsOptimizerState
    : public core::PeriodicOptimizerState {
 private:
  typedef core::PeriodicOptimizerState P;
 private:
  // the model on which the simulation is run and to which all particles are
  // assumed to belong
  base::Pointer<Model> m_;

  int updates_;

  // the types of particles involved in the interaction (type of group I and II)
  // TODO: a bit ugly and ungeneral, we might have mixed types in principle
  InteractionType interaction_type_;

  // maintains a list of nearby particle pairs in a bipartite graph
  IMP::base::PointerMember<IMP::container::CloseBipartitePairContainer>
      close_bipartite_pair_container_;

  // avergae number of times all pairs of particles contacted each other
  // per update round
  Float avg_ncontacts_;

  // average fraction of time that each partilce is in contact
  Float avg_pct_bound_particles_I_;   // for particles in group I
  Float avg_pct_bound_particles_II_;  // particles in group II

  // total number of particles in each group
  Int n_particles_I_;
  Int n_particles_II_;

 public:

  /**
     Constructor

     @param[in] m             the model associated with particlesI / II
     @param[in] interaction_type the types of particles involved in the
                                 interaction, assumed to be the types of
                                 particlesI and II
     @param[in] particlesI    particles from one side of the interaction
     @param[in] particlesII    particles from other side of the interaction
     @param[in] contact_range keep track of particle pairs within that range
     @param[in] slack         slack for updating close particles in appropriate
                              CloseBiparyiyrPairContainer, this affects only
                              performance - touch only if you know why
  */
  BipartitePairsStatisticsOptimizerState(
      Model* m,
      InteractionType interaction_type,  // TODO: remove this from class?
                                         //       a bit ugly and ungeneral
      const ParticlesTemp& particlesI, const ParticlesTemp& particlesII,
      double contact_range = 1.0, double slack = 1.0);

  /**
     returns the particle types of the first and second group of particles,
     respectively. This might change in the future, since types can be mixed
  */
  InteractionType get_interaction_type() const { return interaction_type_; }

  /**
     returns the average number of all pairs of particles that
     touch each other per update
  */
  Float get_average_number_of_contacts() const { return avg_ncontacts_; }

  /**
     returns the average fraction of particles from group I
     that are bound in each update round
  */
  Float get_average_percentage_bound_particles_1() const {
    return avg_pct_bound_particles_I_;
  }

  /**
     returns the average fraction of particles from group II
     that are bound in each update round
  */
  Float get_average_percentage_bound_particles_2() const {
    return avg_pct_bound_particles_II_;
  }

  /**
     return the total number of particles in the first group
  */
  Int get_number_of_particles_1() { return n_particles_I_; }

  /**
     return the total number of particles in the second group
  */
  Int get_number_of_particles_2() { return n_particles_II_; }

  void reset();
  virtual void do_update(unsigned int call_num) IMP_OVERRIDE;
  IMP_OBJECT_METHODS(BipartitePairsStatisticsOptimizerState);
};
IMP_OBJECTS(BipartitePairsStatisticsOptimizerState,
            BipartitePairsStatisticsOptimizerStates);

IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_BIPARTITE_PAIRS_STATISTICS_OPTIMIZER_STATE_H */
