/**
 *  \file initialize_positions.h
 *  \brief description
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_INITIALIZE_POSITIONS_H
#define IMPNPCTRANSPORT_INITIALIZE_POSITIONS_H

#include "npctransport_config.h"
#include "SimulationData.h"

IMPNPCTRANSPORT_BEGIN_NAMESPACE
IMPNPCTRANSPORTEXPORT

/**
   pre-optimize the positions of all diffusing particles in <sd> whose
   coordinates are optimizable, using only chain restraints, excluded
   volumes and bounding box restraints (but not interaction
   restraints).

   The initialization data is dumped to the RMF file `sd->get_rmf_file_name()`
   using dump interval `sd->get_rmf_dump_interval_frames() * 100`, or every
   frame in case that `debug` is true.

   @param sd the simulationd data object containing diffusing particles
   @param extra_restraints a list of additional ad-hoc restraints that will be
                           used only throughout initialization
   @param debug if true, the initialization will dump much more output (e.g.
                every frame to RMF file)
 */
void initialize_positions(SimulationData *sd,
                          const RestraintsTemp &extra_restraints =
                              RestraintsTemp(),
                          bool debug = false);

IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_INITIALIZE_POSITIONS_H */
