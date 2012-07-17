/**
 *  \file simulation_data.cpp
 *  \brief description.
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/automatic_parameters.h>
#include <IMP/atom/estimates.h>
#include <IMP/base/exception.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE
#ifdef IMP_NPC_GOOGLE
using namespace proto2;
#else
using namespace ::google::protobuf;
#endif
#define UPDATE_MIN(name, path)                  \
  min_##name= std::min(min_##name, path().value())
#define UPDATE_MAX(name, path)                  \
  max_##name= std::max(max_##name, path().value())



double get_close_pairs_range(double max_range, double max_range_factor) {
  return max_range*max_range_factor;
}

double get_close_pairs_range(const ::npctransport_proto::Assignment& config) {
  double max_range= config.interaction_range().value();
  UPDATE_MAX(range, config.nonspecific_range);
  double max_range_factor= 1;
  for ( int i=0; i< config.fgs_size(); ++i) {
    UPDATE_MAX(range_factor, config.fgs(i).interaction_range_factor);
  }
  for ( int i=0; i< config.floaters_size(); ++i) {
    UPDATE_MAX(range_factor, config.floaters(i).interaction_range_factor);
  }
  for ( int i=0; i< config.interactions_size(); ++i) {
    if (config.interactions(i).has_interaction_range()) {
      UPDATE_MAX(range, config.interactions(i).interaction_range);
    }
  }
  return get_close_pairs_range(max_range, max_range_factor);
}

double
get_time_step(double max_d_factor,
              double max_k, double min_radius,
              double max_trans_relative_to_radius,
              double time_step_factor) {
  double D=max_d_factor*atom::get_einstein_diffusion_coefficient(min_radius);
  double max_length= max_trans_relative_to_radius * min_radius;
  // binary search between minimal and maximal time steps
  // till they converge near time step that obtains maximal translation
  double ts_max= 1e12, ts_min=0;
  do {
    double mid= .5*(ts_max+ts_min);
    double length= atom::get_diffusion_length(D, max_k, mid);
    if (length > max_length) {
      ts_max=mid;
    } else {
      ts_min=mid;
    }
  } while ((ts_max-ts_min) > .1*ts_max);
  return time_step_factor*ts_min;
}

double get_time_step(const ::npctransport_proto::Assignment& config,
                     double max_trans_relative_to_radius) {
  double time_step_factor= config.time_step_factor().value();
  double max_d_factor=1;
  double min_radius=std::numeric_limits<double>::max();
  double max_k=0;
  UPDATE_MAX(k, config.interaction_k);
  UPDATE_MAX(k, config.backbone_k); // TODO: is this valid for harmonic k?
  UPDATE_MAX(k, config.nonspecific_k);
  UPDATE_MAX(k, config.excluded_volume_k);
  double max_k_factor=1;
  for ( int i=0; i< config.fgs_size(); ++i) {
    UPDATE_MAX(d_factor, config.fgs(i).d_factor);
    UPDATE_MIN(radius, config.fgs(i).radius);
    UPDATE_MAX(k_factor, config.fgs(i).interaction_k_factor);
  }
  for ( int i=0; i< config.floaters_size(); ++i) {
    UPDATE_MAX(d_factor, config.floaters(i).d_factor);
    UPDATE_MIN(radius, config.floaters(i).radius);
    UPDATE_MAX(k_factor, config.floaters(i).interaction_k_factor);
  }
  for ( int i=0; i< config.interactions_size(); ++i) {
    if (config.interactions(i).has_interaction_k()) {
      UPDATE_MAX(k_factor, config.interactions(i).interaction_k);
    }
  }

  return get_time_step(max_d_factor,
                       max_k*max_k_factor,
                       min_radius,
                       max_trans_relative_to_radius,
                       time_step_factor);
}

int get_number_of_frames
(const ::npctransport_proto::Assignment& config, double time_step)
{
  double simulation_time_ns = config.simulation_time_nanosec();
  const double fs_in_ns = 1000000;
  int nframes = std::ceil(simulation_time_ns * fs_in_ns / time_step);
  if(nframes > config.maximal_number_of_frames()){
    IMP_THROW("number of frames " << nframes
              << ", which is required for simulation time" << simulation_time_ns
              << " exceeds the specified maximal value "
              << config.maximal_number_of_frames() ,
              IMP::base::ValueException);
  }
  return nframes;
}

int get_dump_interval_in_frames
(const ::npctransport_proto::Assignment& config, double time_step)
{
  std::cout << "here dump auto" << std::endl;
  const double fs_in_ns = 1000000;
  double ret =
    config.is_dump_interval_in_ns()
    ? config.dump_interval() *  fs_in_ns / time_step
    : config.dump_interval();
  std::cout << "dump interval = " << std::ceil(ret) << " frames, originally "
            << config.dump_interval() << ", time step " << time_step
            << std::endl;;
  return std::ceil(ret);
}

int get_statistics_interval_in_frames
(const ::npctransport_proto::Assignment& config, double time_step)
{
  std::cout << "here stats auto" << std::endl;
  const double fs_in_ns = 1000000;
  double ret =
    config.is_statistics_interval_in_ns()
    ? config.statistics_interval() * fs_in_ns / time_step
    : config.statistics_interval();
  std::cout << "stats interval = " << std::ceil(ret) << " frames, originally "
            << config.statistics_interval() << ", time step " << time_step
            << std::endl;;
  return std::ceil(ret);

}


IMPNPCTRANSPORT_END_NAMESPACE
