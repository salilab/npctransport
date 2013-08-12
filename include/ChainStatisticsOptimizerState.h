/**
 *  \file npctransport/ChainStatisticsOptimizerState.h
 *  \brief description
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_CHAIN_STATISTICS_OPTIMIZER_STATE_H
#define IMPNPCTRANSPORT_CHAIN_STATISTICS_OPTIMIZER_STATE_H

#include "npctransport_config.h"
#include <IMP/Particle.h>
#include <IMP/OptimizerState.h>
#include <IMP/optimizer_state_macros.h>
#include <IMP/core/PeriodicOptimizerState.h>
#include <IMP/npctransport/typedefs.h>
#include <deque>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

/** Compute various statistics of a chain.*/
class IMPNPCTRANSPORTEXPORT ChainStatisticsOptimizerState
    : public core::PeriodicOptimizerState {
 private:
  typedef core::PeriodicOptimizerState P;

  // particles in the chain:
  ParticlesTemp ps_;

  // time series of the positions of particles in the chain:
  std::deque<algebra::Vector3Ds> positions_;

  double get_dt() const;

 public:
  /**
     @param ps the particles being wrapped
     @param periodicity frame interval for statistics, equiv. to set_period(1)
   */
  ChainStatisticsOptimizerState(const ParticlesTemp &ps,
                                unsigned int periodicity = 1);

  double get_correlation_time() const;

  Floats get_diffusion_coefficients() const;

  double get_diffusion_coefficient() const;

  /**
     Resets all the statistics about that chain
  */
  void reset();
  virtual void do_update(unsigned int call_num) IMP_OVERRIDE;
  IMP_OBJECT_METHODS(ChainStatisticsOptimizerState);
};
IMP_OBJECTS(ChainStatisticsOptimizerState, ChainStatisticsOptimizerStates);

IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_CHAIN_STATISTICS_OPTIMIZER_STATE_H */
