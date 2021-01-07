/**
 *  \file PoreRadiusSingletonScore.cpp
 *  \brief XXXX.
 *
 *  Copyright 2007-2021 IMP Inventors. All rights reserved.
 *
 */

#include "IMP/npctransport/PoreRadiusSingletonScore.h"

IMPNPCTRANSPORT_BEGIN_NAMESPACE

ModelObjectsTemp PoreRadiusSingletonScore::do_get_inputs(
    Model *m, const ParticleIndexes &pis) const {
  return IMP::get_particles(m, pis);
}

IMPNPCTRANSPORT_END_NAMESPACE
