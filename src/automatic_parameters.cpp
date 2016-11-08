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

double get_close_pairs_range(const ::npctransport_proto::Assignment& a) {
  double max_range = a.interaction_range().value();
  UPDATE_MAX(range, a.nonspecific_range);
  for (int i = 0; i < a.interactions_size(); ++i) {
    if (a.interactions(i).has_interaction_range()) {
      double k=a.interactions(i).interaction_k().value();
      double range=a.interactions(i).interaction_range().value();
      // TODO: can support is_on, though it's rare
      if(range>0.0 && k>0.0) {
        // compute skewed range if needed
	bool is_orientational=false;
	if(a.interactions(i).has_range_sigma0_deg() &&
	   a.interactions(i).has_range_sigma1_deg()) {
	  if(a.interactions(i).range_sigma0_deg().value()!=0.0 &&
	     a.interactions(i).range_sigma1_deg().value()!=0.0) {
	    is_orientational=true;
	  }
	}
	if(is_orientational){
	  const double pi = 3.1415926535897;
	  double range_sigma0_rad=a.interactions(i).range_sigma0_deg().value()*pi/180.0;
	  double range_sigma1_rad=a.interactions(i).range_sigma1_deg().value()*pi/180.0;
	  std::string type0=a.interactions(i).type0();
	  std::string type1=a.interactions(i).type1();
	  double R0=0.0;
	  double R1=0.0;
	  for(int ii=0; ii<a.floaters_size(); ii++){
	    if(a.floaters(ii).type()==type0){
	      R0=a.floaters(ii).radius().value();
	    }
	    if(a.floaters(ii).type()==type1){
	      R1=a.floaters(ii).radius().value();
	    }
	  }
	  for(int ii=0; ii<a.fgs_size(); ii++){
	    if(a.fgs(ii).type()==type0){
	      R0=a.fgs(ii).radius().value();
	    }
	    if(a.fgs(ii).type()==type1){
	      R1=a.fgs(ii).radius().value();
	    }
	  }
	  double chord0=2*R0*std::sin(range_sigma0_rad/2.0);
	  double chord1=2*R1*std::sin(range_sigma1_rad/2.0);
	  double chord=std::max(chord0,chord1);
	  range=std::max(range, chord);
        }
        // udpate max range
        max_range=std::max(max_range, range);
      }
    }
  }
  double max_range_factor = 0.0001;
  for (int i = 0; i < a.fgs_size(); ++i) {
    UPDATE_MAX(range_factor, a.fgs(i).interaction_range_factor);
  }
  for (int i = 0; i < a.floaters_size(); ++i) {
    UPDATE_MAX(range_factor, a.floaters(i).interaction_range_factor);
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
  double min_radius = std::numeric_limits<double>::max(); // in A
  double min_range = std::numeric_limits<double>::max(); // in A
  double max_k = 0.0; // in kcal/mol/A

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
  }
  for (int i = 0; i < a.floaters_size(); ++i) {
    UPDATE_MAX(d_factor, a.floaters(i).d_factor);
    UPDATE_MIN(radius, a.floaters(i).radius);
  }
  double base_k = a.interaction_k().value();
  double base_range = a.interaction_range().value();

  // go over all interaction and compute their maximal k (in kcal/mol/A)
  // and minimum range (in A)
  for (int i = 0; i < a.interactions_size(); ++i) {
    double k=base_k;
    double range=base_range;
    if (a.interactions(i).has_interaction_k()) {
      k=a.interactions(i).interaction_k().value();
    }
    if (a.interactions(i).has_interaction_range()) {
      range=a.interactions(i).interaction_range().value();
    }
    // factor k and range + retrieve particles radii
    std::string type0=a.interactions(i).type0();
    std::string type1=a.interactions(i).type1();
    double R0(-1.0);
    double R1(-1.0);
    for(int ii=0; ii<a.floaters_size(); ii++){
      if(a.floaters(ii).type()==type0){
        if(a.floaters(ii).interactions().value()==0){
          range=0.0; // skip
          continue;
        }
        k*=a.floaters(ii).interaction_k_factor().value();
        range*=a.floaters(ii).interaction_range_factor().value();
        R0=a.floaters(ii).radius().value();
      }
      if(a.floaters(ii).type()==type1){
        if(a.floaters(ii).interactions().value()==0){
          range=0.0; // skip
          continue;
        }
        k*=a.floaters(ii).interaction_k_factor().value();
        range*=a.floaters(ii).interaction_range_factor().value();
        R1=a.floaters(ii).radius().value();
      }
    }
    for(int ii=0; ii<a.fgs_size(); ii++){
      if(a.fgs(ii).type()==type0){
        if(a.fgs(ii).interactions().value()==0){
          range=0.0; // skip
          continue;
        }
        k*=a.fgs(ii).interaction_k_factor().value();
        range*=a.fgs(ii).interaction_range_factor().value();
        R0=a.fgs(ii).radius().value();
      }
      if(a.fgs(ii).type()==type1){
        if(a.fgs(ii).interactions().value()==0){
          range=0.0; // skip
          continue;
        }
        k*=a.fgs(ii).interaction_k_factor().value();
        range*=a.fgs(ii).interaction_range_factor().value();
        R1=a.fgs(ii).radius().value();
      }
    }
    // compute skewed range if needed
    bool is_orientational=false;
    if(a.interactions(i).has_range_sigma0_deg() &&
       a.interactions(i).has_range_sigma1_deg()) {
      if(a.interactions(i).range_sigma0_deg().value()!=0.0 &&
         a.interactions(i).range_sigma1_deg().value()!=0.0) {
        is_orientational=true;
      }
    }
    if(is_orientational && range > 0.0){
      IMP_USAGE_CHECK(R0>0.0 && R1>0.0,
                      "R0 or R1 could not be found for type0 or type1");
      k*=0.5*range; // the maximal k for this interction
      const double pi = 3.1415926535897;
      double range_sigma0_rad=a.interactions(i).range_sigma0_deg().value()*pi/180.0;
      double range_sigma1_rad=a.interactions(i).range_sigma1_deg().value()*pi/180.0;
      double chord0=2*R0*std::sin(range_sigma0_rad/2.0);
      double chord1=2*R1*std::sin(range_sigma1_rad/2.0);
      double min_chord=std::min(chord0,chord1);
      range=std::min(range, min_chord);
    } // if is_orientation
    if(range>0.0 && k>0.0) {
      max_k=std::max(max_k, k);
      std::cout << "interaction " << i
                << " update range from "
                << min_range;
      min_range=std::min(min_range, range);
      std::cout << " to " << min_range << std::endl;
    }
  } // for interactions(i)

  std::cout << "get_time_step(): "
            << " max_d_factor " << max_d_factor
            << " max-k " << max_k
            << " min-range " << min_range
            << " max-trans-relative-to-R " << max_trans_relative_to_radius
            << " time-step-factor " << time_step_factor
            << std::endl;
  double dT_fs= get_time_step(max_d_factor, max_k, min_radius, min_range,
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
  std::cout << "dump interval = " << ret << " frames, "
            << a.dump_interval_ns() << " ns, time step " << time_step
            << std::endl;
  return ret;
}

int get_statistics_interval_in_frames
( const ::npctransport_proto::Assignment& a,
  double time_step,
  double default_value_ns)
{
  double ns(default_value_ns);
  if(a.has_statistics_interval_ns()) {
    ns = a.statistics_interval_ns() ;
  }
  int ret = get_frames_from_ns(ns, time_step);
  std::cout << "stats interval = "
            << ret  << " frames, converted from "
            << ns << " ns, time step " << time_step
            << std::endl;
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
  std::cout << "output stats interval = "
            << ret  << " frames, converted from "
            << ns << " ns, time step " << time_step
            << std::endl;
  return ret;
}

IMPNPCTRANSPORT_END_NAMESPACE
