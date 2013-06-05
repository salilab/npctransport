/**
 *  \file io.cpp
 *  \brief description.
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
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
#include <IMP/compatibility/set.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE
BodyStatisticsOptimizerState
::BodyStatisticsOptimizerState(Particle*p):
  core::PeriodicOptimizerState("BodyStatisticsOptimizerState%1%"),
  p_(p){
}
void BodyStatisticsOptimizerState::reset() {
  positions_.clear();
  core::PeriodicOptimizerState::reset();
}
double BodyStatisticsOptimizerState::get_dt() const {
  return dynamic_cast<atom::Simulator*>(get_optimizer())
      ->get_maximum_time_step();
}

double BodyStatisticsOptimizerState::
get_correlation_time()const {
  double sum=0;
  int n=0;
  Floats angles;
  for (unsigned int i=0; i< positions_.size(); ++i) {
    double last=0;
    for (unsigned int j=i+1; j < positions_.size(); ++j) {
      algebra::Rotation3D rel= positions_[j].get_rotation()
          /positions_[i].get_rotation();
      double angle=algebra::get_axis_and_angle(rel).second;
      if (i==0) angles.push_back(angle);
      if (angle >1) {
        sum+=get_period()*get_dt()*(j-i-1 + (angle-last));
        ++n;
        break;
      }
      last=angle;
    }
  }
  /*std::cout << n << " events from " << angles
    << " with " << positions_.size() << " samples " << std::endl;*/
  if (n==0) {
    return std::numeric_limits<double>::infinity();
  }
  return sum/n;
}
double BodyStatisticsOptimizerState::get_diffusion_coefficient() const {
  if (positions_.empty()) return 0;
  algebra::Vector3Ds
    displacements(positions_.size()-1);
  for (unsigned int i=1; i< positions_.size(); ++i) {
    displacements[i-1]= positions_[i].get_translation()
        -positions_[i-1].get_translation();
  }
  return atom::get_diffusion_coefficient(displacements,
                                         get_period()*get_dt());
}


void BodyStatisticsOptimizerState
::do_update(unsigned int) {
  positions_.push_back(core::RigidBody(p_).get_reference_frame().
                          get_transformation_to());
  while (positions_.size() > 1000) {
    positions_.pop_front();
  }
}

IMPNPCTRANSPORT_END_NAMESPACE
