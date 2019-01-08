/**
 *  \file io.cpp
 *  \brief description.
 *
 *  Copyright 2007-2019 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/GlobalStatisticsOptimizerState.h>
#include <IMP/npctransport/Statistics.h>
#include <IMP/npctransport/SimulationData.h>
#include <limits>

IMPNPCTRANSPORT_BEGIN_NAMESPACE
GlobalStatisticsOptimizerState::GlobalStatisticsOptimizerState
( IMP::npctransport::Statistics* statistics_manager,
  unsigned int periodicity)
  : P(statistics_manager->get_model(), "BodyStatisticsOptimizerState%1%"),
    statistics_manager_(statistics_manager)
{
  set_period(periodicity);
  reset();
}

void GlobalStatisticsOptimizerState::reset() {
  P::reset();
  mean_energy_= std::numeric_limits<double>::max();
  n_= 0;
}

void GlobalStatisticsOptimizerState::do_update(unsigned int call_num) {
  IMP_UNUSED(call_num);
  double energy=
    statistics_manager_->get_sd()->get_bd()
    ->get_scoring_function()->evaluate(false);
  mean_energy_=
    (energy + mean_energy_*n_)/(n_+1);
  //  IMP_LOG(PROGRESS, "global stats energy=" << energy
  //        " mean energy " << mean_energy_ <<
  //        << " n "  << n_ << std::endl);
  n_++;
}

IMPNPCTRANSPORT_END_NAMESPACE
