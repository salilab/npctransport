/**
 *  \file main.h
 *  \brief Helper functions for executable .cpp files
 *
 *  Copyright 2007-2010 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_INTERNAL_MAIN_H
#define IMPNPCTRANSPORT_INTERNAL_MAIN_H

#include "../npctransport_config.h"
#include "../SimulationData.h"
#include <string>

IMPNPCTRANSPORT_BEGIN_INTERNAL_NAMESPACE

/**
   optimize a system in sd

   @param sd a constructed SimulationData object to be optimized
   @param init_restraints additional restraints to be used ad-hoc
                          during position initialization only
                          (relevant only if init_rmf is "")
   @param quick whether to perform a quick optimization (TODO: this
                right now as sd received 'quick' in its constructor -
                isn't this redundant? shouldn't quick mode be saved
                in sd? or only externally?)
   @param is_initial_optimization whether to do initial optimization of positions
   @param is_equilibration whether to do equilibration BD simulation
   @param is_full_run whetehr to run the full BD simulation
   @param final_config name of a file for saving the final config from each trial
   @param debug if true, do verbose debug outout
*/ //TODO: replace all million flag params with a unified "flags" data structure
IMPNPCTRANSPORTEXPORT void do_main_loop(SimulationData *sd,
                                        const RestraintsTemp &init_restraints,
                                        bool quick,
                                        bool is_initial_optimization,
                                        bool is_equilibration,
                                        bool is_full_run,
                                        std::string final_config,
                                        bool debug_init);
IMPNPCTRANSPORT_END_INTERNAL_NAMESPACE

#endif /* IMPNPCTRANSPORT_INTERNAL_MAIN_H */
