/**
 *  \file ChainStatisticsOptimizerState.cpp
 *  \brief description.
 *
 *  Copyright 2007-2022 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/ChainStatisticsOptimizerState.h>
#include <IMP/npctransport/enums.h>
#include <IMP/npctransport/util.h>
#include <IMP/algebra/geometric_alignment.h>
#include <IMP/atom/estimates.h>
#include <IMP/atom/distance.h>
#include <IMP/atom/Simulator.h>
#include <IMP/core/XYZ.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

ChainStatisticsOptimizerState::ChainStatisticsOptimizerState
( const ParticlesTemp& ps,
  unsigned int periodicity )
  :
  P(ps[0]->get_model(), "ChainStatisticsOptimizerState%1%"),
  ps_(ps),
  mean_rgyr_(-1.0),
  mean_rgyr2_(-1.0),
  mean_end_to_end_(-1.0),
  mean_end_to_end2_(-1.0),
  mean_bond_distance_(-1.0),
  mean_bond_distance2_(-1.0),
  n_(0)
{
  IMP_OBJECT_LOG;
  set_period(periodicity);
  reset();
}

void ChainStatisticsOptimizerState::reset()
{
  IMP_OBJECT_LOG;
  P::reset();
  positions_.clear();
  times_fs_.clear();
  mean_rgyr_= -1.0;
  mean_rgyr2_= -1.0;
  mean_end_to_end_= -1.0;
  mean_end_to_end2_= -1.0;
  mean_bond_distance_= -1.0;
  mean_bond_distance2_= -1.0;
  n_= 0; // resets mean statistics
  core::PeriodicOptimizerState::reset();
}

double ChainStatisticsOptimizerState::get_dt() const {
  return dynamic_cast<atom::Simulator*>(get_optimizer())
      ->get_maximum_time_step();
}

double ChainStatisticsOptimizerState::get_correlation_time() const
{
  IMP_OBJECT_LOG;
  bool IS_DISABLED=true; // DISABLE FOR NOW CAUSE TIME CONSUMING
  if(IS_DISABLED){
    return std::numeric_limits<double>::infinity();
  }
  double sum= 0;
  int n= 0;
  for (unsigned int i= 0; i < positions_.size(); ++i) {
    double last= 0;
    for (unsigned int j= i + 1; j < positions_.size(); ++j) {
      algebra::Transformation3D tr =
          algebra::get_transformation_aligning_first_to_second(positions_[i],
                                                               positions_[j]);
      algebra::Rotation3D rel = tr.get_rotation();
      double angle = algebra::get_axis_and_angle(rel).second;
      if (angle > 1) {
        sum += get_period() * (j - i - 1 + (angle - last)) * get_dt();
        ++n;
        break;
      }
      last = angle;
    }
  }
  IMP_LOG(VERBOSE, n << " correlation events" << std::endl);

  if (n == 0) {
    return std::numeric_limits<double>::infinity();
  }
  return sum / n;
}

Floats
ChainStatisticsOptimizerState::get_local_diffusion_coefficients() const
{
  IMP_OBJECT_LOG;
  if (positions_.empty()) return Floats();
  Vector<algebra::Vector3Ds> displacements(
      positions_[0].size(), algebra::Vector3Ds(positions_.size() - 1));
  for (unsigned int i= 1; i < positions_.size(); ++i) {
    algebra::Transformation3D rel =
        algebra::get_transformation_aligning_first_to_second(positions_[i - 1],
                                                             positions_[i]);
    for (unsigned int j= 0; j < positions_[i].size(); ++j) {
      displacements[j][i - 1] =
          rel.get_transformed(positions_[i - 1][j]) - positions_[i][j];
    }
  }
  Floats ret;
  for (unsigned int i = 0; i < displacements.size(); ++i) {
    ret.push_back(atom::get_diffusion_coefficient(displacements[i],
                                                  get_period() * get_dt()));
  }
  return ret;
}

double ChainStatisticsOptimizerState::get_diffusion_coefficient() const
{
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
  // Compute chain center-of-mass displacement and dT vectors:
  algebra::Vector3Ds positions(positions_.size());
  for (unsigned int i= 0; i < positions_.size(); ++i) {
    positions[i] = std::accumulate(positions_[i].begin(), positions_[i].end(),
                                   algebra::get_zero_vector_d<3>()) /
                   positions_[i].size();
  }
  algebra::Vector3Ds displacements;
  displacements.reserve(n-1);
  std::transform(positions.begin()+1, positions.end(),
                 positions.begin(),
                 std::back_inserter(displacements),
                 std::minus<algebra::Vector3D>());
  IMP::Floats dts;
  dts.reserve(n-1);
  std::transform(times_fs_.begin()+1, times_fs_.end(),
                   times_fs_.begin(),
                 std::back_inserter(dts),
                 std::minus<double>());
  return atom::get_diffusion_coefficient
    (displacements, dts); //  get_period() * get_dt());
}

void ChainStatisticsOptimizerState::do_update(unsigned int) {
  IMP_OBJECT_LOG;
  atom::Simulator* simulator =
    dynamic_cast< atom::Simulator* >( get_optimizer() );
  IMP_USAGE_CHECK( simulator, "Optimizer must be a simulator in order to use "
                   "ChainStatisticsOptimizerState, for time stats" );
  // Get positions and times from current round:
  algebra::Vector3Ds vs;
  for (unsigned int i= 0; i < ps_.size(); ++i) {
    vs.push_back(core::XYZ(ps_[i]).get_coordinates());
  }
  double cur_time_fs = simulator->get_current_time();
  times_fs_.push_back( cur_time_fs );
  positions_.push_back(vs);
  while (positions_.size() > 1000) {
    times_fs_.pop_front();
    positions_.pop_front();
  }
  // Radius of gyration and end-to-end distance of chain/bond:
  double w= 1.0/(++n_);
  double rgyr= atom::get_radius_of_gyration(ps_, false);
  double rgyr2= rgyr*rgyr;
  mean_rgyr_=  w*rgyr  + (1-w)*mean_rgyr_;
  mean_rgyr2_= w*rgyr2 + (1-w)*mean_rgyr2_;
  double end_to_end= core::get_distance( core::XYZ(ps_.front()),
                                         core::XYZ(ps_.back() ) );
  double end_to_end2 = end_to_end*end_to_end;
  mean_end_to_end_=  w*end_to_end +  (1-w)*mean_end_to_end_;
  mean_end_to_end2_= w*end_to_end2 + (1-w)*mean_end_to_end2_;

  double cur_mean_bond_distance(0.0);
  double cur_mean_bond_distance2(0.0);
  for(unsigned int i=1; i<ps_.size(); i++){
    double distance_i= core::get_distance( core::XYZ(ps_[i-1]),
                                           core::XYZ(ps_[i]) );
    cur_mean_bond_distance+=  distance_i/(ps_.size()-1);
    cur_mean_bond_distance2+= (distance_i*distance_i)/(ps_.size()-1);
  }
  mean_bond_distance_= w*cur_mean_bond_distance + (1-w)*mean_bond_distance_;
  mean_bond_distance2_= w*cur_mean_bond_distance2 + (1-w)*mean_bond_distance2_;
}



IMPNPCTRANSPORT_END_NAMESPACE
