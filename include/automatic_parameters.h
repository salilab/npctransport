/**
 *  \file automatic_parameters.h
 *  \brief description
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_AUTOMATIC_PARAMETERS_H
#define IMPNPCTRANSPORT_AUTOMATIC_PARAMETERS_H

#include "npctransport_config.h"
#ifndef SWIG
namespace npctransport_proto {
class Assignment;
}
#endif

IMPNPCTRANSPORT_BEGIN_NAMESPACE
IMPNPCTRANSPORTEXPORT
double get_close_pairs_range(double max_range, double max_range_factor);

/** returns an upper bound on the contact range between any two particles in the
    system */  // TODO: maybe different ranges for different pair types?
IMPNPCTRANSPORTEXPORT
double get_close_pairs_range(const ::npctransport_proto::Assignment& config);

/** Computes the time step size that is required for a stable
    simulation.
    Formally, first computes the time step size that restricts the
    estimated translation size of any particle, at any given simulation step
    to [max_trans_relative_to_radius * min_radius], given the simulation
    parameters max_d_factor and max_k. This base time step is then multiplied by
    time_step_factor.

    @param max_d_factor the maximal diffusion factor of any particle
                        in the system, which factors the Einstein diffusion
                        coefficient)
    @param max_k        the maximal force applied on any particle
                        in the simulation
    @param min_radius   the minimal radius of any particle in the system
    @param max_trans_relative_to_radius the maximal estimated translation
                                        allowed for any particle as fraction of
                                        min_radius (before factoring by
                                        time_step_factor)
    @param time_step_factor multiply final time step in this factor

    @return time step size in femtoseconds required for stable simulation
*/
IMPNPCTRANSPORTEXPORT  // TODO: is max_k also correct for springs?
    double
        get_time_step(double max_d_factor, double max_k, double min_radius,
                      double max_trans_relative_to_radius = 0.1,
                      double time_step_factor = 1.0);

/** computes the time step size that is required for a stable simulation,
    where the translation of any particle at any simulation time step is
    restricted, based on the simulation parameters in config.

    @param config the simulation parameters used to compute the time step size
    @param max_trans_relative_to_radius the maximal estimated translation
                                        allowed for any particle as fraction of
                                        its radius (before factoring by
                                        time step factor, specified in config)

    @return time step size in femtoseconds required for stable simulation
*/
IMPNPCTRANSPORTEXPORT
double get_time_step(const ::npctransport_proto::Assignment& config,
                     double max_trans_relative_to_radius = 0.1);

/**
   computes the number of frames that approximate a certain number of nanoseconds
   using the specified time_step per step in fs

   @param time in ns
   @param time_step step in fs per frame
*/
IMPNPCTRANSPORTEXPORT
int get_frames_from_ns(double ns, double time_step);

/**
   computes the number of frames needed to acheive simulation time
   required in the confiugration time, with time step computed from
   config.

   @param config the simulation parameters
   @param time_step the time step in femtoseconds

   @throw ValueException if maximum number of frames specified in
        config is exceeded */
IMPNPCTRANSPORTEXPORT
int get_number_of_frames(const ::npctransport_proto::Assignment& config,
                         double time_step);

/**
   computes the number of frames for the specified dump_interval,
   which is normalized by the time step if the interval is specified in ns

   @param config the simulation parameters
   @param time_step the time step in femtoseconds
*/
IMPNPCTRANSPORTEXPORT
int get_dump_interval_in_frames(const ::npctransport_proto::Assignment& config,
                                double time_step);

/**
   computes the number of frames for the specified statistics_interval,
   normalized by the time step if the interval is specified in ns

   @param assign the simulation assignment parameters
   @param time_step the time step in femtoseconds
   @param default_value_ns interval in ns to use if none specified in
                           assignment parameters
*/
IMPNPCTRANSPORTEXPORT
int get_statistics_interval_in_frames
( const ::npctransport_proto::Assignment& assign,
  double time_step,
  double default_value_ns = 0.1);

/**
   computes the maximal interval for outputting statistics to output file
   (and to dump order params) in frames, for the output_statistics_interval_ns
   param in assignment, based on specified time_step

   @param assign the simulation assignment parameters
   @param time_step the time step in femtoseconds
   @param default_value_ns interval in ns to use if none specified in
                           assignment parameters
*/
IMPNPCTRANSPORTEXPORT
int get_output_statistics_interval_in_frames
( const ::npctransport_proto::Assignment& assign,
  double time_step,
  double default_value_ns = 1.0);

IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_AUTOMATIC_PARAMETERS_H */
