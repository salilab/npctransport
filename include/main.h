/**
 *  \file main.h
 *  \brief Helper functions for executable .cpp files
 *
 *  Copyright 2007-2021 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_MAIN_H
#define IMPNPCTRANSPORT_MAIN_H

#include "npctransport_config.h"
#include "SimulationData.h"
#include <IMP/base_types.h>
#include <string>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

/**
   initialize and return a simulation data object based on
   program command line parameters

   @throw IMP::IOException if there was any IO problem
*/
IMPNPCTRANSPORTEXPORT
IMP::npctransport::SimulationData *startup(int argc, char *argv[]);

//! remove nup42 and its anchors (also from obstacles)
void remove_Nup42(SimulationData* sd);

//! inflate floater of specified type to new_radius
void inflate_floater
(SimulationData* sd, const std::string floater_name, const float new_radius);

//! change box size sd to specified box size and update output file
void reset_box_size(SimulationData* sd, double box_size);

/** Run simulation using preconstructed SimulationData object sd.

    @param sd SimulationData object to optimize
    @param init_restraints ad-hoc restraints during initialization only
*/
IMPNPCTRANSPORTEXPORT
void do_main_loop(SimulationData *sd, const RestraintsTemp &init_restraints);

IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_MAIN_H */
