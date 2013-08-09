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
#include <IMP/npctransport/typedefs.h>
#include <IMP/atom/Simulator.h>
#include <deque>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

/**
     Maintains transport statistics about a particle p in a z-axis
     aligned channel,
 */
class IMPNPCTRANSPORTEXPORT ParticleTransportStatisticsOptimizerState
    : public core::PeriodicOptimizerState {
 private:
  typedef core::PeriodicOptimizerState P;
  Particle *p_;             // the particle
  Float bottom_z_, top_z_;  // channel boundaries on z-axis
  base::WeakPointer<IMP::atom::Simulator> owner_;
  unsigned int n_transports_up_;        // from bottom to top of channel
  unsigned int n_transports_down_;      // from top to bottom of channOBel
  Floats transport_time_points_in_ns_;  // simulation time points of each
                                        // transport event
  // MOVED NEXT LINE TO DECORATOR (TODO: perhaps also other stats)
  // bool is_last_entry_from_top_; // last time p entered channel
  //                               // was top or bottom
  bool is_reset_;

  Particle *get_particle() const { return p_; }

 public:
  /**
     Initiates transport statistics about a particle p in a z-axis aligned
     channel,
     whose bottom and top are at z-coordinates bottom_z and top_z, respectively.

     @param p the particle, assumed to be decorated as a RigidBody
     @param bottom_z the z coordinate of the channel bottom
     @param top_z the z coordinate of the channel top
     @param owner a simulator that is moving this particle and can provide it
     with
                  time information, or nullptr
   */
  ParticleTransportStatisticsOptimizerState(
      Particle *p, Float bottom_z, Float top_z,
      base::WeakPointer<IMP::atom::Simulator> owner = nullptr);

  //! sets a simulator that moves this particle and can provide simulation time
  //! information about it, or IMP::nullptr if none
  void set_owner(base::WeakPointer<IMP::atom::Simulator> owner) {
    owner_ = owner;
  }

  //! returns the simulator that was declared in the constructor or by
  //set_owner()
  //! to moves this particle, and provide simulation time information about it.
  base::WeakPointer<IMP::atom::Simulator> get_owner() const { return owner_; }

  /**
      Returns the number of times the particle crossed the channel
      from its bottom to its top
  */
  unsigned int get_n_transports_up() const { return n_transports_up_; }

  /** Returns the number of times the particle crossed the channel
      from its top to its bottom
  */
  unsigned int get_n_transports_down() const { return n_transports_down_; }

  /** Returns the number of times the particle crossed the channel
      from any one side to the other
  */
  unsigned int get_total_n_transports() const {
    return n_transports_up_ + n_transports_down_;
  }

  /**
     returns a list of all the simulation time points in nanoseconds
     when a transport even has occured (according to the owner of this
     OptimizerState). An empty list is returned if no
     transport events are known or if there was no owner when they occured
   */
  Floats const &get_transport_time_points_in_ns() const {
    return transport_time_points_in_ns_;
  }

  /** resets the number of transports statistics to 0 */
  void reset();
  virtual void do_update(unsigned int call_num) IMP_OVERRIDE;
  IMP_OBJECT_METHODS(ParticleTransportStatisticsOptimizerState);
};

IMP_OBJECTS(ParticleTransportStatisticsOptimizerState,
            ParticleTransportStatisticsOptimizerStates);

IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_PARTICLE_TRANSPORT_STATISTICS_OPTIMIZER_STATE_H */
