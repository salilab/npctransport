/**
 *  \file simulation_data.cpp
 *  \brief description.
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/automatic_parameters.h>
#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 6)
#pragma GCC diagnostic push
#endif
#pragma GCC diagnostic ignored "-Wsign-compare"
#include "npctransport.pb.h"
#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 6)
#pragma GCC diagnostic pop
#endif

#include <IMP/atom/estimates.h>
#include <IMP/exception.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE
using namespace IMP_NPCTRANSPORT_PROTOBUF_NAMESPACE;
#define UPDATE_MIN(name, path) min_##name = std::min(min_##name, path().value())
#define UPDATE_MAX(name, path) max_##name = std::max(max_##name, path().value())

double get_close_pairs_range(double max_range, double max_range_factor) {
  // squared cause range factor is applied once for each interacting partner
  return max_range * max_range_factor * max_range_factor;
}

double get_close_pairs_range(const ::npctransport_proto::Assignment& config) {
  double max_range = config.interaction_range().value();
  UPDATE_MAX(range, config.nonspecific_range);
  for (int i = 0; i < config.interactions_size(); ++i) {
    if (config.interactions(i).has_interaction_range()) {
      UPDATE_MAX(range, config.interactions(i).interaction_range);
    }
  }
  double max_range_factor = 0.0001;
  for (int i = 0; i < config.fgs_size(); ++i) {
    UPDATE_MAX(range_factor, config.fgs(i).interaction_range_factor);
  }
  for (int i = 0; i < config.floaters_size(); ++i) {
    UPDATE_MAX(range_factor, config.floaters(i).interaction_range_factor);
  } // TODO: add obstacles?!
  return get_close_pairs_range(max_range, max_range_factor);
}

double get_time_step(double max_d_factor, double max_k, double min_radius,
                     double max_trans_relative_to_radius,
                     double time_step_factor) {
  double D =
      max_d_factor * atom::get_einstein_diffusion_coefficient(min_radius);
  double max_length = max_trans_relative_to_radius * min_radius;
  // binary search between minimal and maximal time steps
  // till they converge near time step that obtains maximal translation
  double ts_max = 1e12, ts_min = 0;
  do {
    double mid = .5 * (ts_max + ts_min);
    double length = atom::get_diffusion_length(D, max_k, mid);
    if (length > max_length) {
      ts_max = mid;
    } else {
      ts_min = mid;
    }
  } while ((ts_max - ts_min) > .1 * ts_max);
  return time_step_factor * ts_min;
}

double get_time_step(const ::npctransport_proto::Assignment& a,
                     double max_trans_relative_to_radius)
{
  // NOTE: this is not a tight bound - a tight bound would involve
  // computing time step for each interactions (factored by relevant
  // constants and accounting for particle radii and interaction range
  // / k). So this is a conservative time step choice in that respect.
  double time_step_factor = a.time_step_factor().value();
  double max_d_factor = 1.0;
  double min_radius = std::numeric_limits<double>::max();
  double max_k = 0.0;
  double max_k_factor = 1.0;
  double min_radius_factor = 1.0;
  UPDATE_MAX(k, a.interaction_k);
  UPDATE_MAX(k, a.backbone_k);  // TODO: is this valid for harmonic k?
  UPDATE_MAX(k, a.nonspecific_k);
  UPDATE_MIN(radius, a.nonspecific_range)
  UPDATE_MAX(k, a.excluded_volume_k);
  for (int i = 0; i < a.fgs_size(); ++i) {
    UPDATE_MAX(d_factor, a.fgs(i).d_factor);
    UPDATE_MIN(radius, a.fgs(i).radius);
    UPDATE_MAX(k_factor, a.fgs(i).interaction_k_factor);
    UPDATE_MIN(range_factor, a.fgs(i).interaction_range_factor);
  }
  for (int i = 0; i < a.floaters_size(); ++i) {
    UPDATE_MAX(d_factor, a.floaters(i).d_factor);
    UPDATE_MIN(radius, a.floaters(i).radius);
    UPDATE_MAX(k_factor, a.floaters(i).interaction_k_factor);
    UPDATE_MIN(range_factor, a.floaters(i).interaction_range_factor);
  }
  for (int i = 0; i < a.interactions_size(); ++i) {
    if (a.interactions(i).has_interaction_k() && a.interactions(i).has_interaction_range()) {
      if(a.interactions(i).interaction_k().value()>0.0 &&
         a.interactions(i).interaction_range().value()>0.0) {
        UPDATE_MAX(k_factor, a.interactions(i).interaction_k);
        std::cout << "interaction #" << i
                  << " range: " << a.interactions(i).interaction_range().value()
                  << std::endl;
        UPDATE_MIN(radius, a.interactions(i).interaction_range);
      }
    }
  }

  std::cout << "get_time_step(): "
            << " max_d_factor " << max_d_factor
            << " max-k " << max_k
            << " max-k-factor " << max_k_factor << " min_radius " << min_radius
            << " min-range/radius-factor " << min_radius_factor
            << " max-trans-relative-to-R " << max_trans_relative_to_radius
            << " time-step-factor " << time_step_factor
            << std::endl;
  double dT_fs= get_time_step(max_d_factor, max_k * max_k_factor, min_radius * min_radius_factor,
                       max_trans_relative_to_radius, time_step_factor);
  std::cout << "dT = " << dT_fs << " [fs]" << std::endl;
  return dT_fs;

}

int get_frames_from_ns(double ns, double time_step) {
  const double fs_in_ns = 1000000;
  int frames = std::ceil(ns * fs_in_ns / time_step);
  if (frames <= 0) frames = 1;  // make sure at least every frame
  return frames;
}


int get_number_of_frames(const ::npctransport_proto::Assignment& a,
                         double time_step) {
  int ret = get_frames_from_ns(a.simulation_time_ns(), time_step);
  if (ret > a.maximal_number_of_frames()) {
    IMP_THROW(
        "number of frames " << ret << ", which is required for simulation time"
                            << a.simulation_time_ns()
                            << " exceeds the specified maximal value "
                            << a.maximal_number_of_frames(),
        IMP::ValueException);
  }
  return ret;
}


int get_dump_interval_in_frames(const ::npctransport_proto::Assignment& a,
                                double time_step) {
  int ret = get_frames_from_ns(a.dump_interval_ns(), time_step);
  IMP_LOG(VERBOSE, "dump interval = " << ret << " frames, "
            << a.dump_interval_ns() << " ns, time step " << time_step
            << std::endl);
  ;
  return ret;
}

int get_statistics_interval_in_frames
( const ::npctransport_proto::Assignment& a,
  double time_step,
  double default_value_ns)
{
  double ns = default_value_ns;
  if(a.has_statistics_interval_ns()) {
    ns = a.statistics_interval_ns() ;
  }
  int ret = get_frames_from_ns(ns, time_step);
  IMP_LOG(PROGRESS, "stats interval = "
          << ret  << " frames, converted from "
            << ns << " ns, time step " << time_step
            << std::endl);
  return ret;
}

int get_output_statistics_interval_in_frames
( const ::npctransport_proto::Assignment& a,
  double time_step,
  double default_value_ns)
{
  double ns = default_value_ns;
  if(a.has_output_statistics_interval_ns()) {
    ns = a.output_statistics_interval_ns() ;
  }
  int ret = get_frames_from_ns(ns, time_step);
  IMP_LOG(PROGRESS, "output stats interval = "
          << ret  << " frames, converted from "
          << ns << " ns, time step " << time_step
          << std::endl
          );
  return ret;
}

IMPNPCTRANSPORT_END_NAMESPACE
