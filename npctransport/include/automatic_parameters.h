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

IMPNPCTRANSPORTEXPORT
double get_close_pairs_range(const ::npctransport::Assignment& config);

IMPNPCTRANSPORTEXPORT
double get_time_step(double time_step_factor, double max_d_factor,
              double max_k, double min_radius,
                     double min_range) ;
IMPNPCTRANSPORTEXPORT
double get_time_step(const ::npctransport::Assignment& config);

IMPNPCTRANSPORT_END_NAMESPACE

#endif  /* IMPNPCTRANSPORT_AUTOMATIC_PARAMETERS_H */
