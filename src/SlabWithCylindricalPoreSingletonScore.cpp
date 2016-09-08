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
(double thickness, double radius, double k)
    : thickness_(thickness),
      radius_(radius),
      k_(k),
      top_(thickness / 2.0),
      bottom_(-thickness / 2.0),
      midZ_(0.0)
{
    IMP_LOG_PROGRESS("Constructing a slab with cylindrical pore singleton score"
                     << " from z="  << bottom_ << " to z=" << top_
                     << "; R=" << radius_ << "; k=" << k << std::endl);

}

ModelObjectsTemp SlabWithCylindricalPoreSingletonScore::do_get_inputs(
    Model *m, const ParticleIndexes &pis) const {
  return IMP::get_particles(m, pis);
}

IMPNPCTRANSPORT_END_NAMESPACE
