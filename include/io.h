/**
 *  \file npctransport/io.h
 *  \brief description
 *
 *  Copyright 2007-2021 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_IO_H
#define IMPNPCTRANSPORT_IO_H

#include "npctransport_config.h"
#include <IMP/macros.h>
#include <RMF/FileHandle.h>
#include <IMP/display/Writer.h>
#include <IMP/OptimizerState.h>
#include <IMP/algebra/VectorD.h>
#include <RMF/FileHandle.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

class SimulationData;
IMPNPCTRANSPORTEXPORT void write_geometry(const ParticlesTemp &kaps,
                                          const algebra::Sphere3Ds &kap_sites,
                                          const ParticlesTemp &chains,
                                          const algebra::Sphere3Ds &chain_sites,
                                          const RestraintsTemp &rs,
                                          display::Writer *out);

IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_IO_H */
