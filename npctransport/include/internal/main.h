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
   @param links
   @param quick whether to perform a quick optimization (TODO: this
                right now as sd received 'quick' in its constructor -
                isn't this redundant? shouldn't quick mode be saved
                in sd? or only externally?)
   @param final_config name of a file for saving the final config
   @param debug_init if true, do the initialization in a verbose debug
                     mode, with lots of output
   @param init_rmf if of length one or more, initialize the coordinates
                   of the system from the last frame of the specified RMF file,
                   likely that from a previous run using the same run params
 */
IMPNPCTRANSPORTEXPORT void do_main_loop(SimulationData *sd,
                                        const ParticlePairsTemp &links,
                                        bool quick,
                                        std::string final_config,
                                        bool debug_init,
                                        std::string init_rmf="");
IMPNPCTRANSPORT_END_INTERNAL_NAMESPACE

#endif /* IMPNPCTRANSPORT_INTERNAL_MAIN_H */
