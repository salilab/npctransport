/**
 *  \file automatic_parameters.h
 *  \brief description
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_AUTOMATIC_PARAMETERS_H
#define IMPNPCTRANSPORT_AUTOMATIC_PARAMETERS_H

#include "npctransport_config.h"
#ifdef IMP_NPC_GOOGLE
#include "third_party/npc/module/data/npctransport.pb.h"
#else
#include "npctransport.pb.h"
#endif
IMPNPCTRANSPORT_BEGIN_NAMESPACE
IMPNPCTRANSPORTEXPORT
double get_close_pairs_range(double max_range, double max_range_factor);

/** returns an upper bound on the contact range between any two particles in the
    system */ // TODO: maybe different ranges for different pair types?
IMPNPCTRANSPORTEXPORT
double get_close_pairs_range(const ::npctransport::Assignment& config);

/** Computes the time step size that is required for a stable
    simulation (=restricting the size of each particle motion at each step),
    given certain simulation parameters.

    @param time_step_factor multiply final time step in this factor
    @param max_d_factor the maximal diffusion facotr of any particle
                        in the system
    @param max_k        the maximal force applied on any particle
                        in the simulation
    @param min_radius   the minimal radius of any particle in the system
    @param min_range    the minimal attraction / repulsion range between
                        any two particles in the simulation
    @param max_trans_relative_to_radius required maximal translation of
                                         any particle (per time_step_factor)
                                         as fraction of min_radius,
    @param max_trans_relative_to_range required maximal translation of
                                         any particle (per time_step_factor)
                                         as fraction of min_range

    @return time step size in femtoseconds required for stable simulation
*/
IMPNPCTRANSPORTEXPORT // TODO: is max_k also correct for springs?
double get_time_step(double time_step_factor, double max_d_factor,
                     double max_k, double min_radius,
                     double min_range,
                     double max_trans_relative_to_radius = 0.1,
                     double max_trans_relative_to_range = 0.3) {
) ;
IMPNPCTRANSPORTEXPORT
double get_time_step(const ::npctransport::Assignment& config);

IMPNPCTRANSPORT_END_NAMESPACE

#endif  /* IMPNPCTRANSPORT_AUTOMATIC_PARAMETERS_H */
