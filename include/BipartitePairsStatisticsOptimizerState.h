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
//#include <IMP/optimizer_state_macros.h>
#include <IMP/core/PeriodicOptimizerState.h>
#include <IMP/container/CloseBipartitePairContainer.h>
#include <IMP/npctransport/typedefs.h>
#include <boost/unordered_set.hpp>
#include <deque>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

class Statistics;

/** Track the interaction between pairs from one group of particles
    with particles from another group, within some specified contat
    range
*/
class IMPNPCTRANSPORTEXPORT BipartitePairsStatisticsOptimizerState
    : public core::PeriodicOptimizerState {
 private:
  typedef core::PeriodicOptimizerState P;
  typedef IMP_KERNEL_LARGE_ORDERED_SET<ParticleIndex> t_particle_index_ordered_set; // got to be ordered
  typedef IMP_KERNEL_LARGE_ORDERED_SET<ParticleIndexPair> t_particle_index_pair_ordered_set;
 private:
  WeakPointer<IMP::npctransport::Statistics> statistics_manager_;

  bool is_reset_;

  int n_updates_;
  double time_ns_; // simulation time in ns
  double stats_time_ns_; // time statistics were gathered for miscs
  double off_stats_time_ns_; // time when off-rate calculations were gathered,
                            // requiring presence of bounds, weighted by #contacts
  double off_I_stats_time_ns_; // weighted by number of bounds I
  double off_II_stats_time_ns_; // weighted by number of bounds II
  double on_stats_time_ns_; // time when on-rate calcualtions were gathered,
                            // requiring presence of unbounds
  double on_I_stats_time_ns_; // time when on-rate calcualtions were gathered
                              // for particles I, weighted by their number
  double on_II_stats_time_ns_; // time when on-rate calcualtions were gathered
                              // for particles II, weighted by their number

  // the types of particles involved in the interaction (type of group I and II)
  // TODO: a bit ugly and ungeneral, we might have mixed types in principle
  InteractionType interaction_type_;

  // maintains a list of nearby particle pairs in a bipartite graph
  IMP::PointerMember<IMP::container::CloseBipartitePairContainer>
      close_bipartite_pair_container_;

  // range considered as contact (without slack of CloseBipartitePairContainer)
  double range_;

  // avergae number of times all pairs of particles contacted each other
  // per update round
  double avg_ncontacts_;

  // list of bound particles of each type + list of their contacts after
  // last round of update
  t_particle_index_ordered_set bounds_I_;
  t_particle_index_ordered_set bounds_II_;
  t_particle_index_pair_ordered_set contacts_;

  // Average since last reset:
  double avg_pct_bound_particles_I_; // particles in group I
  double avg_pct_bound_particles_II_;  // particles in group II
  double avg_off_per_contact_per_ns_;
  double avg_off_per_bound_I_per_ns_;
  double avg_off_per_bound_II_per_ns_;
  double avg_on_per_unbound_I_per_ns_;
  double avg_on_per_unbound_II_per_ns_;
  double avg_on_per_missing_contact_per_ns_;

  // total number of all (bound) particles in each group after last udpate()
  unsigned int n_particles_I_;
  unsigned int n_particles_II_;
  unsigned int n_bounds_I_;
  unsigned int n_bounds_II_;

  // theoretical number of possible (unordered) contacts between particles
  unsigned int n_possible_contacts_;

 public:

  /**
     Constructor

     @param[in] m             the model associated with particlesI / II
     @param[in] interaction_type the pair of interacting particle types, assumed to
                                 be the types of all particles in particlesI and
                                 particlesII, respectively.
     @param[in] particlesI    particles from one side of the interaction
     @param[in] particlesII    particles from other side of the interaction
     @param[in] contact_range define a contact for particle pairs whose sphere distances
                              are within contact_range in [A]
     @param[in] slack         slack for updating close particles in appropriate
                              CloseBiparyiyrPairContainer, this affects only
                              performance - touch only if you know what you're
                              doing
  */
  BipartitePairsStatisticsOptimizerState
    (WeakPointer<IMP::npctransport::Statistics> statistics_manager,
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
     returns the average number of interacting pairs of particles
     per update
  */
  double get_average_number_of_contacts() const { return avg_ncontacts_; }

  /** returns the average off rate per bound complex per nanosecond
      since last reset()
  */
  double get_average_off_per_contact_per_ns() const
  { return avg_off_per_contact_per_ns_; }

  /** returns the average off rate per bound particles of type I
      since last reset()
  */
  double get_average_off_per_bound_I_per_ns() const
  { return avg_off_per_bound_I_per_ns_; }

  /** returns the average off rate per bound particles of type II
      since last reset()
  */
  double get_average_off_per_bound_II_per_ns() const
  { return avg_off_per_bound_II_per_ns_; }

  /** returns the average on rate per missing possible contact
      (out of all theoretically possible contacts between particlesI
      and particles II)
  */
  double get_average_on_per_missing_contact_per_ns() const
  { return avg_on_per_missing_contact_per_ns_; }

  /** returns the average on rate per unbound particles of type I
      since last resetI()
  */
  double get_average_on_per_unbound_I_per_ns() const
  { return avg_on_per_unbound_I_per_ns_; }

  /** returns the average on rate per unbound particles of type II
      since last resetII()
  */
  double get_average_on_per_unbound_II_per_ns() const
  { return avg_on_per_unbound_II_per_ns_; }

  /**
     returns the average fraction of particles from group I
     that are bound in each update round
  */
  double get_average_percentage_bound_particles_1() const {
    return avg_pct_bound_particles_I_;
  }

  /**
     returns the average fraction of particles from group II
     that are bound in each update round
  */
  double get_average_percentage_bound_particles_2() const {
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


  /** restart accumulation of all averages in the next time
      that update() is called
  */
  void reset();

  double get_misc_stats_period_ns() const
  { return stats_time_ns_; };

  double get_off_stats_period_ns() const
  { return off_stats_time_ns_; };

  double get_off_I_stats_period_ns() const
  { return off_I_stats_time_ns_; };

    double get_off_II_stats_period_ns() const
  { return off_II_stats_time_ns_; };

    double get_on_stats_period_ns() const
  { return on_stats_time_ns_; };

  double get_on_I_stats_period_ns() const
  { return on_I_stats_time_ns_; };

  double get_on_II_stats_period_ns() const
  { return on_II_stats_time_ns_; };


 protected:
  virtual void do_update(unsigned int call_num) IMP_OVERRIDE;


 private:

  /**
     update fraction bound with n1 bound particles of type I;
     and n2 bound particles of type II, assuming old_updates
     for the old average
   */
  //  void update_fraction_bound(unsigned int n1,
  //                            unsigned int n2,
  //                           unsigned int old_updates);


 public:
  IMP_OBJECT_METHODS(BipartitePairsStatisticsOptimizerState);
};
IMP_OBJECTS(BipartitePairsStatisticsOptimizerState,
            BipartitePairsStatisticsOptimizerStates);

IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_BIPARTITE_PAIRS_STATISTICS_OPTIMIZER_STATE_H */
