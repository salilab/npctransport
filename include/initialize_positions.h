/**
 *  \file initialize_positions.h
 *  \brief description
 *
 *  Copyright 2007-2018 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_INITIALIZE_POSITIONS_H
#define IMPNPCTRANSPORT_INITIALIZE_POSITIONS_H

#include "npctransport_config.h"
#include "SimulationData.h"

IMPNPCTRANSPORT_BEGIN_NAMESPACE
IMPNPCTRANSPORTEXPORT

/**
   pre-optimize the positions of all diffusing particles in 'sd' whose
   coordinates are optimizable, using only chain restraints, excluded
   volumes and bounding box restraints (but not interaction
   restraints).

   The initialization data is dumped to the RMF file `sd->get_rmf_file_name()`
   using dump interval `sd->get_rmf_dump_interval_frames() * 100`, or every
   frame in case that `debug` is true.

   @param sd the SimulationData object containing diffusing particles
   @param extra_restraints a list of additional ad-hoc restraints that will be
                           used only throughout initialization
   @param debug if true, the initialization will dump much more output (e.g.
                every frame to RMF file)
   @param short_init_factor a factor between >0 and 1 for decreasing
                            the number of optimization cycles at each
                            round
   @param is_disable_randomize if true, do not initially randomize particle positions,
                               essentially performing an extended relaxation from
                               the starting coordinates
   @param are_fgs_pre_initialized if true, do not try to pre-optimize FGs before adding diffusers
 */
void initialize_positions
( SimulationData *sd,
  const RestraintsTemp &extra_restraints =
  RestraintsTemp(),
  bool debug = false,
  double short_init_factor = 1.0,
  bool is_disable_randomize = false,
  bool are_fgs_pre_initialized = false );


IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_INITIALIZE_POSITIONS_H */
