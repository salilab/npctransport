/**
 *  \file SlabWithToroidalPorePairScore.cpp
 *  \brief a score for a slab with a toroidal pore
 *
 *  Copyright 2006-15 Sali Lab. All rights reserved.
 *
 */

#include "IMP/npctransport/SlabWithToroidalPorePairScore.h"

IMPNPCTRANSPORT_BEGIN_NAMESPACE

SlabWithToroidalPorePairScore
::SlabWithToroidalPorePairScore
(double k) :
  k_(k)
{
  IMP_LOG_PROGRESS("Constructing a slab with toroidal pore singleton score"
                   << "; k=" << k << std::endl);
}

ModelObjectsTemp
SlabWithToroidalPorePairScore
::do_get_inputs
(Model *m, const ParticleIndexes &pis) const
{
  return IMP::get_particles(m, pis);
}

IMPNPCTRANSPORT_END_NAMESPACE
