/**
 *  \file npctransport/BodyStatisticsOptimizerState.h
 *  \brief description
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_BODY_STATISTICS_OPTIMIZER_STATE_H
#define IMPNPCTRANSPORT_BODY_STATISTICS_OPTIMIZER_STATE_H

#include "npctransport_config.h"
#include <IMP/Particle.h>
#include <IMP/algebra/Transformation3D.h>
#include <IMP/OptimizerState.h>
//#include <IMP/optimizer_state_macros.h>
#include <IMP/core/PeriodicOptimizerState.h>
#include <IMP/npctransport/typedefs.h>
#include <deque>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

/** Track the rotational correlation time of a rigid body particle*/
/** The correlation with at most the last 100 updates is tracked*/
class IMPNPCTRANSPORTEXPORT BodyStatisticsOptimizerState
    : public core::PeriodicOptimizerState {
 private:
  typedef core::PeriodicOptimizerState P;
  Particle *p_;
  std::deque<algebra::Transformation3D> positions_;
  Particle *get_particle() const { return p_; }
  void add_orientation(algebra::Rotation3D rot) { positions_.push_back(rot); }
  double get_dt() const;

 public:
  /**
     @param p the particle being wrapped
     @param periodicity frame interval for statistics, equiv. to set_period(1)
   */
  BodyStatisticsOptimizerState(Particle *p, unsigned int periodicity=1);
  double get_correlation_time() const;
  double get_diffusion_coefficient() const;
  void reset();
  virtual void do_update(unsigned int call_num) IMP_OVERRIDE;
  IMP_OBJECT_METHODS(BodyStatisticsOptimizerState);
};
IMP_OBJECTS(BodyStatisticsOptimizerState, BodyStatisticsOptimizerStates);

IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_BODY_STATISTICS_OPTIMIZER_STATE_H */
