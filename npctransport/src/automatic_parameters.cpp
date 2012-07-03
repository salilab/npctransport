/**
 *  \file simulation_data.cpp
 *  \brief description.
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/automatic_parameters.h>
#include <IMP/atom/estimates.h>

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

double get_close_pairs_range(const ::npctransport::Assignment& config) {
  double max_range= config.interaction_range().value();
  UPDATE_MAX(range, config.nonspecific_range);
  double max_range_factor= 1;
  for ( int i=0; i< config.fgs_size(); ++i) {
    UPDATE_MAX(range_factor, config.fgs(i).interaction_range_factor);
  }
  for ( int i=0; i< config.floaters_size(); ++i) {
    UPDATE_MAX(range_factor, config.floaters(i).interaction_range_factor);
  }

  return get_close_pairs_range(max_range, max_range_factor);
}

double
get_time_step(double time_step_factor, double max_d_factor,
              double max_k, double min_radius,
              double min_range) {
  double D=max_d_factor*atom::get_einstein_diffusion_coefficient(min_radius);
  double scale= std::min(.1*min_radius, .3*min_range);
  // binary search between ts_max and ts_min
  double ts_max= 1e12, ts_min=0;
  do {
    double mid= .5*(ts_max+ts_min);
    double length= atom::get_diffusion_length(D, max_k, mid);
    if (length > scale) {
      ts_max=mid;
    } else {
      ts_min=mid;
    }
  } while ((ts_max-ts_min) > .1*ts_max);
  return time_step_factor*ts_min;
}

double get_time_step(const ::npctransport::Assignment& config) {
  double time_step_factor= config.time_step_factor().value();
  double max_d_factor=1;
  double min_radius=std::numeric_limits<double>::max();
  double min_range= config.interaction_range().value();
  UPDATE_MIN(range, config.nonspecific_range);
  double min_range_factor= std::numeric_limits<double>::max();
  double max_k=0;
  UPDATE_MAX(k, config.interaction_k);
  UPDATE_MAX(k, config.backbone_k);
  UPDATE_MAX(k, config.nonspecific_k);
  UPDATE_MAX(k, config.excluded_volume_k);
  double max_k_factor=1;
  for ( int i=0; i< config.fgs_size(); ++i) {
    UPDATE_MAX(d_factor, config.fgs(i).d_factor);
    UPDATE_MIN(radius, config.fgs(i).radius);
    UPDATE_MIN(range_factor, config.fgs(i).interaction_range_factor);
    UPDATE_MAX(k_factor, config.fgs(i).interaction_k_factor);
  }
  for ( int i=0; i< config.floaters_size(); ++i) {
    UPDATE_MAX(d_factor, config.floaters(i).d_factor);
    UPDATE_MIN(radius, config.floaters(i).radius);
    UPDATE_MIN(range_factor, config.floaters(i).interaction_range_factor);
    UPDATE_MAX(k_factor, config.floaters(i).interaction_k_factor);
  }
  for ( int i=0; i< config.interactions_size(); ++i) {
    if (config.interactions(i).has_interaction_range()) {
      UPDATE_MIN(range_factor, config.interactions(i)
                 .interaction_range);
    }
    if (config.interactions(i).has_interaction_k()) {
      UPDATE_MAX(k_factor, config.interactions(i).interaction_k);
    }
  }



  return get_time_step(time_step_factor, max_d_factor,
                       max_k*max_k_factor,
                       min_radius, min_range);
}

IMPNPCTRANSPORT_END_NAMESPACE
