/**
 *  \file npctransport/ParticleTransportStatisticsOptimizerState.h
 *  \brief description
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_PARTICLE_TRANSPORT_STATISTICS_OPTIMIZER_STATE_H
#define IMPNPCTRANSPORT_PARTICLE_TRANSPORT_STATISTICS_OPTIMIZER_STATE_H

#include "npctransport_config.h"
#include <IMP/Particle.h>
#include <IMP/OptimizerState.h>
#include <IMP/optimizer_state_macros.h>
#include <IMP/core/PeriodicOptimizerState.h>
#include <IMP/core/periodic_optimizer_state_macros.h>
#include <IMP/npctransport/typedefs.h>
#include <deque>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

/**
     Maintains transport statistics about a particle p in a z-axis aligned channel,
 */
class IMPNPCTRANSPORTEXPORT ParticleTransportStatisticsOptimizerState:
public core::PeriodicOptimizerState {
  Particle *p_; // the particle
  Float bottom_z_, top_z_; // channel boundaries on z-axis
  Float prev_z_; // particle z in previous round
  Float cur_z_; // particle z in this round
  Int n_transports_up_; // from bottom to top of channel
  Int n_transports_down_; // from top to bottom of channOBel
  Int n_entries_bottom_; // times particle entered channel from bottom
  Int n_entries_top_; // times particle entered channel from top
  bool is_last_entry_from_top_; // last time p entered channel was top or bottom

  Particle *get_particle() const {return p_;}
 public:
  /**
     Initiates transport statistics about a particle p in a z-axis aligned channel,
     whose bottom and top are at z-coordinates bottom_z and top_z, respectively.

     @param p the particle, assumed to be decorated as a RigidBody
     @param bottom_z the z coordinate of the channel bottom
     @param top_z the z coordinate of the channel top
   */
  ParticleTransportStatisticsOptimizerState(Particle *p,
                                            Float bottom_z,
                                            Float top_z);

  /** Returns the number of times the particle crossed the channel
      from its bottom to its top
  */
  Float get_n_transports_up()
  { return n_transports_up_; }

  /** Returns the number of times the particle crossed the channel
      from its top to its bottom
  */
  Float get_n_transports_down()
  { return n_transports_down_; }

  /** Returns the number of times the particle crossed the channel
      from any one side to the other
  */
  Float get_total_n_transports()
  { return n_transports_up_ + n_transports_down_; }

  /** resets the number of transports statistics to 0 */
  void reset();

  IMP_CORE_PERIODIC_OPTIMIZER_STATE(ParticleTransportStatisticsOptimizerState);
};

IMP_OBJECTS(ParticleTransportStatisticsOptimizerState,
            ParticleTransportStatisticsOptimizerStates);


IMPNPCTRANSPORT_END_NAMESPACE

#endif  /* IMPNPCTRANSPORT_PARTICLE_TRANSPORT_STATISTICS_OPTIMIZER_STATE_H */
