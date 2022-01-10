/**
 *  \file npctransport/GlobalStatisticsOptimizerState.h
 *  \brief description
 *
 *  Copyright 2007-2022 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_GLOBAL_STATISTICS_OPTIMIZER_STATE_H
#define IMPNPCTRANSPORT_GLOBAL_STATISTICS_OPTIMIZER_STATE_H

#include "npctransport_config.h"
#include <IMP/OptimizerState.h>
//#include <IMP/optimizer_state_macros.h>
#include <IMP/core/PeriodicOptimizerState.h>
#include <IMP/npctransport/typedefs.h>
#include <deque>

IMPNPCTRANSPORT_BEGIN_NAMESPACE
//TODO: add tester
class Statistics;

/** Track the rotational correlation time of a rigid body particle*/
/** The correlation with at most the last 100 updates is tracked*/
class IMPNPCTRANSPORTEXPORT GlobalStatisticsOptimizerState
    : public core::PeriodicOptimizerState {
 private:
  typedef core::PeriodicOptimizerState P;
  WeakPointer<IMP::npctransport::Statistics> statistics_manager_;
  double mean_energy_;
  int n_;

 public:
  /**
     @param statistics_manager an optional statistical manager to which statistical updates are sent
                               (of e.g. zr-distribution that are collectively gathered)
     @param periodicity frame interval for statistics, equiv. to set_period(1)
   */
  GlobalStatisticsOptimizerState
    (IMP::npctransport::Statistics* statistics_manager = nullptr,
     unsigned int periodicity=1);

  void reset() IMP_OVERRIDE;

  virtual void do_update(unsigned int call_num) IMP_OVERRIDE;

  double get_mean_energy() {
    return mean_energy_;
  }

  IMP_OBJECT_METHODS(GlobalStatisticsOptimizerState);
};
IMP_OBJECTS(GlobalStatisticsOptimizerState, GlobalStatisticsOptimizerStates);

IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_GLOBAL_STATISTICS_OPTIMIZER_STATE_H */
