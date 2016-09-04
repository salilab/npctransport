/**
 *  \file SlabWithToroidalPoreSingletonScore.cpp
 *  \brief a score for a slab with a toroidal pore
 *
 *  Copyright 2006-15 Sali Lab. All rights reserved.
 *
 */

#include "IMP/npctransport/SlabWithToroidalPoreSingletonScore.h"

IMPNPCTRANSPORT_BEGIN_NAMESPACE

SlabWithToroidalPoreSingletonScore::SlabWithToroidalPoreSingletonScore
(double radius, double bottom, double top, double k) :
  bottom_(bottom),
  top_(top),
  midZ_((top+bottom)/2.0),
  R_(radius),
  r_((top-bottom)/2.0),
  k_(k)
{}

SlabWithToroidalPoreSingletonScore::SlabWithToroidalPoreSingletonScore
    (double radius, double thickness, double k) :
      SlabWithToroidalPoreSingletonScore(radius, -thickness, thickness, k)
{}

ModelObjectsTemp
SlabWithToroidalPoreSingletonScore::do_get_inputs
(Model *m, const ParticleIndexes &pis) const
{
  return IMP::get_particles(m, pis);
}

IMPNPCTRANSPORT_END_NAMESPACE
