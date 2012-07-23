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
void initialize_positions(SimulationData *sd,
                          const ParticlePairsTemp &extra_links
                          = ParticlePairsTemp(),
                          bool debug=false);


IMPNPCTRANSPORT_END_NAMESPACE

#endif  /* IMPNPCTRANSPORT_INITIALIZE_POSITIONS_H */
