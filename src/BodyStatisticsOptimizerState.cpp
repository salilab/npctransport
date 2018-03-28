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
#include <IMP/npctransport/enums.h>

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
  P::reset();
  positions_.clear();
  times_fs_.clear();
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

namespace{
  struct subtract_translations {
    algebra::Vector3D operator()(algebra::Transformation3D const &a,
                                 algebra::Transformation3D const &b)
    {
      return a.get_translation() - b.get_translation();
    }
  };
};

double BodyStatisticsOptimizerState::get_diffusion_coefficient() const {
    IMP_OBJECT_LOG;
  // Checks:
  unsigned int n= positions_.size();
  IMP_USAGE_CHECK(times_fs_.size() == n,
                  "Length of times and positions lists is expected to be equal");
  if (n<2){
    return 0;
  }
  if (times_fs_.front() - times_fs_.back() == 0.0) {
    return 0;
  }
  // Compute displacements and dT vectors:
  algebra::Vector3Ds displacements;
  displacements.reserve(n-1);
  std::transform(positions_.begin()+1, positions_.end(),
                 positions_.begin(),
                 std::back_inserter(displacements),
                 subtract_translations());
  IMP::Floats dts;
  dts.reserve(n-1);
  std::transform(times_fs_.begin()+1, times_fs_.end(),
                   times_fs_.begin(),
                 std::back_inserter(dts),
                 std::minus<double>());
  return atom::get_diffusion_coefficient
    (displacements, dts); //  get_period() * get_dt());
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
  atom::Simulator* simulator =
    dynamic_cast< atom::Simulator* >( get_optimizer() );
  double cur_time_ns = simulator->get_current_time();
  times_fs_.push_back( cur_time_ns );
  positions_.push_back(core::RigidBody(p_)
                           .get_reference_frame().get_transformation_to());
  while (positions_.size() > 1000) {
    times_fs_.pop_front();
    positions_.pop_front();
  }
}

IMPNPCTRANSPORT_END_NAMESPACE
