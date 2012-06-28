/**
 *  \file npctransport/statistics_optimizer_states.h
 *  \brief description
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_STATISTICS_OPTIMIZER_STATES_H
#define IMPNPCTRANSPORT_STATISTICS_OPTIMIZER_STATES_H

#include "npctransport_config.h"
#include <IMP/Particle.h>
#include <IMP/algebra/Transformation3D.h>
#include <IMP/OptimizerState.h>
#include <IMP/optimizer_state_macros.h>
#include <IMP/container/CloseBipartitePairContainer.h>
#include <IMP/npctransport/typedefs.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

/** Track the rotational correlation time of a rigid body*/
class IMPNPCTRANSPORTEXPORT BodyStatisticsOptimizerState:
public OptimizerState {
  Particle *p_;
  algebra::Transformation3Ds positions_;
  Particle *get_particle() const {return p_;}
  void add_orientation(algebra::Rotation3D rot) {
    positions_.push_back(rot);
  }
  double get_dt() const;
 public:
  BodyStatisticsOptimizerState(Particle *p);
  double get_correlation_time()const;
  double get_diffusion_coefficient() const;
  IMP_PERIODIC_OPTIMIZER_STATE(BodyStatisticsOptimizerState);
};
IMP_OBJECTS(BodyStatisticsOptimizerState,
            BodyStatisticsOptimizerStates);

/** Compute various statistics of a chain.*/
class IMPNPCTRANSPORTEXPORT ChainStatisticsOptimizerState:
public OptimizerState {
  ParticlesTemp ps_;
  base::Vector<algebra::Vector3Ds > positions_;
  double get_dt() const;
 public:
  ChainStatisticsOptimizerState(const ParticlesTemp &p);
  double get_correlation_time()const;
  Floats get_diffusion_coefficients() const;
  double get_diffusion_coefficient() const;
  IMP_PERIODIC_OPTIMIZER_STATE(ChainStatisticsOptimizerState);
};
IMP_OBJECTS(ChainStatisticsOptimizerState,
            ChainStatisticsOptimizerStates);

/** Track the interaction between pairs from one group of particles
    with particles from another group, within some specified contat
    range
*/
class IMPNPCTRANSPORTEXPORT BipartitePairsStatisticsOptimizerState:
public OptimizerState {
  // the model on which the simulation is run and to which all particles are
  // assumed to belong
  Pointer<Model> m_;

  // the types of particles involved in the interaction (type of group I and II)
  // TODO: a bit ugly and ungeneral, we might have mixed types in principle
  InteractionType interaction_type_;

  // maintains a list of nearby particle pairs in a bipartite graph
  IMP::OwnerPointer<IMP::container::CloseBipartitePairContainer>
    close_bipartite_pair_container_;

  // avergae number of times all pairs of particles contacted each other
  // per update round
  Float avg_ncontacts_;

  // average fraction of time that each partilce is in contact
  Float avg_pct_bound_particles_I_; // for particles in group I
  Float avg_pct_bound_particles_II_; // particles in group II

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
  BipartitePairsStatisticsOptimizerState
    (Model* m,
     InteractionType interaction_type, // TODO: remove this from class?
                                        //       a bit ugly and ungeneral
     const ParticlesTemp& particlesI,
     const ParticlesTemp& particlesII,
     double contact_range = 1.0,
     double slack = 1.0 );

  // returns the particle types of the first and second group of particles,
  // respectively. This might change in the future, since types can be mixed
  InteractionType get_interaction_type() const
  { return interaction_type_; }

  // returns the average number of all pairs of particles that
  //touch each other per update
  Float get_avg_ncontacts() const
  { return avg_ncontacts_; }

  // returns the average fraction of particles from group I
  // that are bound in each update round
  Float get_avg_pct_bound_particles_I() const
  { return avg_pct_bound_particles_I_; }

  // returns the average fraction of particles from group II
  // that are bound in each update round
  Float get_avg_pct_bound_particles_II() const
  { return avg_pct_bound_particles_II_; }

  // return the total number of particles in the first group
  Int get_n_particles_groupI()
  { return n_particles_I_; }

  // return the total number of particles in the second group
  Int get_n_particles_groupII()
  { return n_particles_II_; }

  IMP_PERIODIC_OPTIMIZER_STATE(BipartitePairsStatisticsOptimizerState);
};
IMP_OBJECTS(BipartitePairsStatisticsOptimizerState,
            BipartitePairsStatisticsOptimizerStates);



IMPNPCTRANSPORT_END_NAMESPACE

#endif  /* IMPNPCTRANSPORT_STATISTICS_OPTIMIZER_STATES_H */
