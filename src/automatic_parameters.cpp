/**
 *  \file simulation_data.cpp
 *  \brief description.
 *
 *  Copyright 2007-2020 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/automatic_parameters.h>
#include <IMP/atom/estimates.h>
#include <IMP/exception.h>
#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 6)
#pragma GCC diagnostic push
#endif
#if defined(__GNUC__)
#pragma GCC diagnostic ignored "-Wsign-compare"
#endif
#include "npctransport.pb.h"
#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 6)
#pragma GCC diagnostic pop
#endif
#include <limits>

IMPNPCTRANSPORT_BEGIN_NAMESPACE
using namespace IMP_NPCTRANSPORT_PROTOBUF_NAMESPACE;
#define UPDATE_MIN(name, path) min_##name = std::min(min_##name, path().value())
#define UPDATE_MAX(name, path) max_##name = std::max(max_##name, path().value())

double get_close_pairs_range(double max_range, double max_range_factor) {
  // squared cause range factor is applied once for each interacting partner
  double ret_value( max_range * max_range_factor * max_range_factor );
  IMP_LOG(VERBOSE, " Close pairs range: " << ret_value << std::endl);
  return ret_value;
}

namespace {
  //! get radius of particle with named type
  //! in assignment
  //! (assuming only one particle type has that name!)
  IMP_UNUSED_FUNCTION
  double get_radius_of_type
  (const ::npctransport_proto::Assignment a,
   std::string type)
  {
    for(int i=0; i<a.floaters_size(); i++){
      if(a.floaters(i).type()==type){
        return a.floaters(i).radius().value();
      }
    }
    for(int i=0; i<a.fgs_size(); i++){
      if(a.fgs(i).type()==type){
        return a.fgs(i).radius().value();
      }
    }
    for(int i=0; i<a.obstacles_size(); i++){
      if(a.obstacles(i).type()==type){
        return a.obstacles(i).radius().value();
      }
    }
    return 0.0;
  }

  //! get range factor of particle with named type
  //! in assignment
  //! (assuming only one particle type has that name!)
  double get_range_factor_of_type
  (const ::npctransport_proto::Assignment a,
   std::string type)
  {
    for(int i=0; i<a.floaters_size(); i++){
      if(a.floaters(i).type()==type){
        return a.floaters(i).interaction_range_factor().value();
      }
    }
    for(int i=0; i<a.fgs_size(); i++){
      if(a.fgs(i).type()==type){
        return a.fgs(i).interaction_range_factor().value();
      }
    }
    for(int i=0; i<a.obstacles_size(); i++){
      if(a.obstacles(i).type()==type){
        return a.obstacles(i).interaction_range_factor().value();
      }
    }
    return 1.0;
  }

  double get_close_pairs_range_for_interaction
  (const ::npctransport_proto::Assignment a, int interaction_id)
  {
    ::npctransport_proto::Assignment::InteractionAssignment interaction
      (a.interactions(interaction_id));
    std::string type0=interaction.type0();
    std::string type1=interaction.type1();
    double range_factor_0= get_range_factor_of_type(a, type0);
    double range_factor_1= get_range_factor_of_type(a, type1);
    if(!interaction.has_interaction_range()){
      double default_range = a.interaction_range().value();
      return default_range*range_factor_0*range_factor_1;
      return 0.0;
    }
    double range=interaction.interaction_range().value();
    {
      bool is_on= interaction.is_on().value()>0;
      double k=interaction.interaction_k().value();
      double epsilon(std::numeric_limits<double>::epsilon());
      if(range < epsilon || k < epsilon || !is_on)
        {
          return 0.0;
        }
    }
    // Handle orientational case if needed (NOTE: THIS IS DISABLED CAUSE ORIENTATIONAL
    // CASE MAY RESULT IN LONGER SITE-SITE INTERACTION RANGE BUT NOT SPHERE-SPHERE RANGE)
    // bool is_orientational=false;
    // if(interaction.has_range_sigma0_deg() &&
    //    interaction.has_range_sigma1_deg()) {
    //   is_orientational=
    //     (interaction.range_sigma0_deg().value()!=0.0 &&
    //      interaction.range_sigma1_deg().value()!=0.0);
    // }
    // if(is_orientational)
    //   {
    //     const double pi = 3.1415926535897;
    //     double range_sigma0_rad=interaction.range_sigma0_deg().value()*pi/180.0;
    //     double range_sigma1_rad=interaction.range_sigma1_deg().value()*pi/180.0;
    //     double R0= get_radius_of_type(a, type0);
    //     double R1= get_radius_of_type(a, type1);
    //     double chord0=2*R0*std::sin(range_sigma0_rad/2.0);
    //     double chord1=2*R1*std::sin(range_sigma1_rad/2.0);
    //     double chord=std::max(chord0,chord1);
    //     range=std::max(range, chord);
    //   }
    // TODO: do obstacles need special treatment?!
    return range*range_factor_0*range_factor_1;
  }

  //! returns true if type matches fg_assignment.type(), with
  //! or without suffixes from fg_assignments's suffix list
  bool is_fg_type_match
  ( std::string type,
    const ::npctransport_proto::Assignment_FGAssignment& fg_assignment)
  {
    std::string fg_type= fg_assignment.type();
    if(type==fg_type){
      return true;
    }
    for(int i=0; i<fg_assignment.type_suffix_list_size(); i++){
      std::string fg_type_i= fg_type+fg_assignment.type_suffix_list(i);
      if(type==fg_type_i){
        return true;
      }
    }
    return false;
  }
}; // namespace {}

double get_close_pairs_range(const ::npctransport_proto::Assignment& a) {
  double max_range = std::numeric_limits<double>::epsilon();
  UPDATE_MAX(range, a.nonspecific_range);
  for (int i = 0; i < a.interactions_size(); ++i) {
    max_range= std::max(max_range,
                        get_close_pairs_range_for_interaction(a, i));
  }
  return max_range;
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
    bool is_found0= false;
    bool is_found1= false;
    for(int ii=0; ii<a.floaters_size(); ii++){
      if(a.floaters(ii).type()==type0){
        if(a.floaters(ii).interactions().value()==0){
          range=0.0; // skip
          continue;
        }
        k*=a.floaters(ii).interaction_k_factor().value();
        range*=a.floaters(ii).interaction_range_factor().value();
        R0=a.floaters(ii).radius().value();
        is_found0= true;
      }
      if(a.floaters(ii).type()==type1){
        if(a.floaters(ii).interactions().value()==0){
          range=0.0; // skip
          continue;
        }
        k*=a.floaters(ii).interaction_k_factor().value();
        range*=a.floaters(ii).interaction_range_factor().value();
        R1=a.floaters(ii).radius().value();
        is_found1= true;
      }
    }
    for(int ii=0; ii<a.fgs_size(); ii++){
      if(is_fg_type_match(type0, a.fgs(ii))){
        if(a.fgs(ii).interactions().value()==0){
          range=0.0; // skip
          continue;
        }
        k*=a.fgs(ii).interaction_k_factor().value();
        range*=a.fgs(ii).interaction_range_factor().value();
        R0=a.fgs(ii).radius().value();
        is_found0= true;
      }
      if(is_fg_type_match(type1, a.fgs(ii))){
        if(a.fgs(ii).interactions().value()==0){
          range=0.0; // skip
          continue;
        }
        k*=a.fgs(ii).interaction_k_factor().value();
        range*=a.fgs(ii).interaction_range_factor().value();
        R1=a.fgs(ii).radius().value();
        is_found1= true;
      }
    }
    if(!is_found0 || !is_found1){
      IMP_LOG(TERSE, "get_time_step() - ignoring interaction "
                << type0 << "-" << type1
                << " since no particles of one these types are defined"
                << std::endl);
      continue;
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
    if(is_orientational && range > 0.0 && is_found0 && is_found1){
      IMP_USAGE_CHECK(R0>0.0 && R1>0.0,
                      "R0 or R1 could not be found for type " << type0
                      << " or type " << type1 << std::endl);
      k*=0.5*range; // the maximal force for this interction
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
      min_range=std::min(min_range, range);
    }
  } // for interactions(i)

  IMP_LOG(VERBOSE, "get_time_step(): "
            << " max_d_factor " << max_d_factor
            << " max-k " << max_k
            << " min-range " << min_range
            << " max-trans-relative-to-R " << max_trans_relative_to_radius
            << " time-step-factor " << time_step_factor
            << std::endl);
  double dT_fs= get_time_step(max_d_factor, max_k, min_radius, min_range,
                              max_trans_relative_to_radius, time_step_factor);
  IMP_LOG(VERBOSE, "dT = " << dT_fs << " [fs]" << std::endl);
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
  IMP_LOG(VERBOSE, "stats interval = "
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
  IMP_LOG(VERBOSE, "output stats interval = "
            << ret  << " frames, converted from "
            << ns << " ns, time step " << time_step
            << std::endl);
  return ret;
}

IMPNPCTRANSPORT_END_NAMESPACE
