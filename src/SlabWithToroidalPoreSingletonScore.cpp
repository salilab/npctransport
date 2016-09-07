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
(double slab_bottom, double slab_top, double radius, double k) :
  bottom_(slab_bottom),
  top_(slab_top),
  midZ_((slab_top+slab_bottom)/2.0),
  R_(radius),
  r_((slab_top-slab_bottom)/2.0),
  k_(k)
{}

SlabWithToroidalPoreSingletonScore::SlabWithToroidalPoreSingletonScore
(double slab_thickness, double radius, double k) :
  SlabWithToroidalPoreSingletonScore
  (-slab_thickness, slab_thickness, radius, k)
{}

ModelObjectsTemp
SlabWithToroidalPoreSingletonScore::do_get_inputs
(Model *m, const ParticleIndexes &pis) const
{
  return IMP::get_particles(m, pis);
}

IMPNPCTRANSPORT_END_NAMESPACE
