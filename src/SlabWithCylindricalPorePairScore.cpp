/**
 *  \file SlabWithCylindricalPoreSingletonScore.cpp
 *  \brief XXXX.
 *
 *  Copyright 2007-8 Sali Lab. All rights reserved.
 *
 */

#include "IMP/npctransport/SlabWithCylindricalPoreSingletonScore.h"

IMPNPCTRANSPORT_BEGIN_NAMESPACE

SlabWithCylindricalPoreSingletonScore::SlabWithCylindricalPoreSingletonScore
(double k)
    : k_(k)
{
    IMP_LOG_PROGRESS("Constructing a slab with cylindrical pore singleton score"
		     << "; k=" << k << std::endl);

}

ModelObjectsTemp SlabWithCylindricalPoreSingletonScore::do_get_inputs(
    Model *m, const ParticleIndexes &pis) const {
  return IMP::get_particles(m, pis);
}

IMPNPCTRANSPORT_END_NAMESPACE
