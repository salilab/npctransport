/**
 *  \file io.cpp
 *  \brief description.
 *
 *  Copyright 2007-2018 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/BodyStatisticsOptimizerState.h>
#include <IMP/core/rigid_bodies.h>
#include <IMP/algebra/geometric_alignment.h>
#include <IMP/atom/estimates.h>
#include <IMP/atom/Simulator.h>
#include <IMP/container/ClosePairContainer.h>
#include <IMP/container/ConsecutivePairContainer.h>
#include <IMP/container/ListSingletonContainer.h>
#include <IMP/npctransport/Statistics.h>
#include <IMP/npctransport/SimulationData.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE
BodyStatisticsOptimizerState::BodyStatisticsOptimizerState
( Particle* p,
  IMP::npctransport::Statistics* statistics_manager,
  unsigned int periodicity)
  : P(p->get_model(), "BodyStatisticsOptimizerState%1%"), p_(p),
        statistics_manager_(statistics_manager)
{
  set_period(periodicity);
}

void BodyStatisticsOptimizerState::reset() {
  positions_.clear();
  core::PeriodicOptimizerState::reset();
}

double BodyStatisticsOptimizerState::get_dt() const {
  return dynamic_cast<atom::Simulator*>(get_optimizer())
      ->get_maximum_time_step();
}

double BodyStatisticsOptimizerState::get_correlation_time() const {
  bool IS_DISABLED=true; // DISABLE FOR NOW CAUSE TIME CONSUMING
  if(IS_DISABLED){
    return std::numeric_limits<double>::infinity();
  }
  double sum = 0;
  int n = 0;
  Floats angles;
  for (unsigned int i = 0; i < positions_.size(); ++i) {
    double last = 0;
    for (unsigned int j = i + 1; j < positions_.size(); ++j) {
      algebra::Rotation3D rel =
          positions_[j].get_rotation() / positions_[i].get_rotation();
      double angle = algebra::get_axis_and_angle(rel).second;
      if (i == 0) angles.push_back(angle);
      if (angle > 1) {
        sum += get_period() * get_dt() * (j - i - 1 + (angle - last));
        ++n;
        break;
      }
      last = angle;
    }
  }
  /*std::cout << n << " events from " << angles
    << " with " << positions_.size() << " samples " << std::endl;*/
  if (n == 0) {
    return std::numeric_limits<double>::infinity();
  }
  return sum / n;
}
double BodyStatisticsOptimizerState::get_diffusion_coefficient() const {
  if (positions_.empty()) return 0;
  algebra::Vector3Ds displacements(positions_.size() - 1);
  for (unsigned int i = 1; i < positions_.size(); ++i) {
    displacements[i - 1] =
        positions_[i].get_translation() - positions_[i - 1].get_translation();
  }
  return atom::get_diffusion_coefficient(displacements,
                                         get_period() * get_dt());
}

void
BodyStatisticsOptimizerState
::update_particle_type_zr_distribution_map() {
  if (statistics_manager_ == nullptr){
    return;
  }
  if(statistics_manager_->get_sd()->get_is_xyz_hist_stats()){
    statistics_manager_->update_particle_type_xyz_distribution_map(p_);
  } else {
    statistics_manager_->update_particle_type_zr_distribution_map(p_);
  }
}


void BodyStatisticsOptimizerState::do_update(unsigned int) {
  this->update_particle_type_zr_distribution_map();
  positions_.push_back(core::RigidBody(p_)
                           .get_reference_frame().get_transformation_to());
  while (positions_.size() > 1000) {
    positions_.pop_front();
  }
}

IMPNPCTRANSPORT_END_NAMESPACE
