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
IMPNPCTRANSPORTEXPORT void do_main_loop(SimulationData *sd,
                                        const ParticlePairsTemp &links,
                                        bool quick,
                                        std::string final_config,
                                        bool debug_init);
IMPNPCTRANSPORT_END_INTERNAL_NAMESPACE

#endif /* IMPNPCTRANSPORT_INTERNAL_MAIN_H */
