/**
 *  \file ChainStatisticsOptimizerState.cpp
 *  \brief description.
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/ChainStatisticsOptimizerState.h>
#include <IMP/algebra/geometric_alignment.h>
#include <IMP/atom/estimates.h>
#include <IMP/atom/Simulator.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

ChainStatisticsOptimizerState
::ChainStatisticsOptimizerState(const ParticlesTemp&p):
  core::PeriodicOptimizerState("ChainStatisticsOptimizerState%1%"),
  ps_(p) {
}
void ChainStatisticsOptimizerState::reset() {
  positions_.clear();
  core::PeriodicOptimizerState::reset();
}
double ChainStatisticsOptimizerState::get_dt() const {
  return dynamic_cast<atom::Simulator*>(get_optimizer())
      ->get_maximum_time_step();
}
double ChainStatisticsOptimizerState::
get_correlation_time()const {
  double sum=0;
  int n=0;
  for (unsigned int i=0; i< positions_.size(); ++i) {
    double last=0;
    for (unsigned int j=i+1; j < positions_.size(); ++j) {
      algebra::Transformation3D tr
          =algebra::get_transformation_aligning_first_to_second(positions_[i],
                                                                positions_[j]);
      algebra::Rotation3D rel = tr.get_rotation();
      double angle=algebra::get_axis_and_angle(rel).second;
      if (angle >1) {
        sum+=get_period()*(j-i-1 + (angle-last)) *get_dt();
        ++n;
        break;
      }
      last=angle;
    }
  }
  std::cout << n << " events" << std::endl;

  if (n==0) {
    return std::numeric_limits<double>::infinity();
  }
  return sum/n;
}

Floats ChainStatisticsOptimizerState::get_diffusion_coefficients() const {
  if (positions_.empty()) return Floats();
  base::Vector<algebra::Vector3Ds >
    displacements(positions_[0].size(),
                  algebra::Vector3Ds( positions_.size()-1));
  for (unsigned int i=1; i< positions_.size(); ++i) {
    algebra::Transformation3D rel
        = algebra::get_transformation_aligning_first_to_second(positions_[i-1],
                                                               positions_[i]);
    for (unsigned int j=0; j < positions_[i].size(); ++j) {
      displacements[j][i-1]= rel.get_transformed(positions_[i-1][j])
        - positions_[i][j];
    }
  }
  Floats ret;
  for (unsigned int i=0; i < displacements.size(); ++i) {
    ret.push_back(atom::get_diffusion_coefficient(displacements[i],
                                                  get_period()*get_dt()));
  }
  return ret;
}

void ChainStatisticsOptimizerState
::do_update(unsigned int) {
  algebra::Vector3Ds vs;
  for (unsigned int i=0; i< ps_.size(); ++i) {
    vs.push_back(core::XYZ(ps_[i]).get_coordinates());
  }
  positions_.push_back(vs);
  while (positions_.size() > 1000) {
    positions_.pop_front();
  }
}

double ChainStatisticsOptimizerState::get_diffusion_coefficient() const {
  if (positions_.empty()) return 0;
  algebra::Vector3Ds positions(positions_.size());
  for (unsigned int i=0; i< positions_.size(); ++i) {
    positions[i]= std::accumulate(positions_[i].begin(), positions_[i].end(),
                           algebra::get_zero_vector_d<3>())
        /positions_[i].size();
  }
  algebra::Vector3Ds
    displacements(positions.size()-1);
  for (unsigned int i=1; i< positions_.size(); ++i) {
    displacements[i-1]= positions[i] -positions[i-1];
  }
  return atom::get_diffusion_coefficient(displacements,
                                         get_period()*get_dt());
}


IMPNPCTRANSPORT_END_NAMESPACE
