/**
 *  \file SeparateSingletonModifier.cpp
 *  \brief XXXX.
 *
 *  Copyright 2007-8 Sali Lab. All rights reserved.
 *
 */

#include "IMP/npctransport/SlabSingletonScore.h"

IMPNPCTRANSPORT_BEGIN_NAMESPACE

SlabSingletonScore::SlabSingletonScore(double thickness, double radius,
                                       double k)
    : thickness_(thickness),
      radius_(radius),
      k_(k),
      top_(thickness / 2.0),
      bottom_(-thickness / 2.0),
      midZ_(0.0) {}

ModelObjectsTemp SlabSingletonScore::do_get_inputs(
    Model *m, const ParticleIndexes &pis) const {
  return IMP::get_particles(m, pis);
}

IMPNPCTRANSPORT_END_NAMESPACE
