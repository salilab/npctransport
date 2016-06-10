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

double get_time_step(double max_d_factor, double max_k, double min_radius, double min_range,
                     double max_trans_relative_to_radius,
                     double time_step_factor) {
  double D =
      max_d_factor * atom::get_einstein_diffusion_coefficient(min_radius);
  double max_length = max_trans_relative_to_radius * std::min(min_radius, min_range);
  // binary search between minimal and maximal time steps
  // till they converge near time step that obtains maximal translation
  double ts_max = 1e12, ts_min = 0;
  do {
    double ts_mid = .5 * (ts_max + ts_min);
    // maximal estimated movement due to either force or random diffusion
    double length = std::max(atom::get_diffusion_length(D, max_k, ts_mid),
                             atom::get_diffusion_length(D, ts_mid));
    if (length > max_length) {
      ts_max = ts_mid;
    } else {
      ts_min = ts_mid;
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
  double min_range = std::numeric_limits<double>::max();
  double max_k = 0.0;
  double max_k_factor = 1.0;
  double min_range_factor = 1.0;
  UPDATE_MAX(k, a.interaction_k);
  UPDATE_MAX(k, a.backbone_k);  // TODO: is this valid for harmonic k?
  if(a.nonspecific_range().value()>0.0 &&
     a.nonspecific_k().value()>0.0) {
    UPDATE_MAX(k, a.nonspecific_k);
    UPDATE_MIN(range, a.nonspecific_range);
  }
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
      double k=a.interactions(i).interaction_k().value();
      double range=a.interactions(i).interaction_range().value();
      if(range>0.0 && k>0.0) {
	bool is_skewed=false;
	if(a.interactions(i).has_k_tangent_skew() &&
	   a.interactions(i).has_range_tangent_skew()) {
	  if(a.interactions(i).k_tangent_skew().value()>0.0 &&
	     a.interactions(i).range_tangent_skew().value()>0.0) {
	    is_skewed=true;
	  }
	}
	if(is_skewed){
	  // normalize for skew and compute maximal force k (k*range/2.0) in skewed interaction force field
	  double k_tangent_skew=a.interactions(i).k_tangent_skew().value();
	  double range_tangent_skew=a.interactions(i).range_tangent_skew().value();
	  double range1=range*std::sqrt(range_tangent_skew);
	  double range2=range/std::sqrt(range_tangent_skew);
	  double k1=k*std::sqrt(k_tangent_skew)*range1/2.0;
	  double k2=k/std::sqrt(k_tangent_skew)*range2/2.0;
	  k=    std::max(k1,k2);
	  range=std::min(range1,range2);
	  std::cout << "Skewed interaction detected - k1*range1/2.0 " << k1 << " k2*range2/2.0 " << k2 
		    << " range1 " << range1 << " range 2 " << range2 << std::endl;
	}
	max_k=    std::max(max_k,k);
	min_range=std::min(min_range,range);
	
      }
    }
  }
  
  std::cout << "get_time_step(): "
            << " max_d_factor " << max_d_factor
            << " max-k " << max_k
            << " max-k-factor " << max_k_factor << " min_radius " << min_radius
            << " min-range " << min_range
            << " min-range-factor " << min_range_factor
            << " max-trans-relative-to-R " << max_trans_relative_to_radius
            << " time-step-factor " << time_step_factor
            << std::endl;
  double dT_fs= get_time_step(max_d_factor, max_k * max_k_factor, min_radius, min_range * min_range_factor * min_range_factor,
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
