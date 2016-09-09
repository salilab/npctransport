/**
 *  \file SlabWithToroidalPoreSingletonScore.cpp
 *  \brief a score for a slab with a toroidal pore
 *
 *  Copyright 2006-15 Sali Lab. All rights reserved.
 *
 */

#include "IMP/npctransport/SlabWithToroidalPoreSingletonScore.h"

IMPNPCTRANSPORT_BEGIN_NAMESPACE

SlabWithToroidalPoreSingletonScore
::SlabWithToroidalPoreSingletonScore
(double slab_thickness, double radius, double k, double horizontal_minor_radius)
  :
  midZ_(0.0),
  R_(radius),
  rv_(0.5*slab_thickness),
  rh_(horizontal_minor_radius),
  k_(k),
  bottom_(midZ_-rv_),
  top_(midZ_+rv_)
{}

SlabWithToroidalPoreSingletonScore
::SlabWithToroidalPoreSingletonScore
(double slab_thickness, double radius, double k) :
  midZ_(0.0),
  R_(radius),
  rv_(0.5*slab_thickness),
  rh_(rv_),
  k_(k),
  bottom_(midZ_-rv_),
  top_(midZ_+rv_)
{}

ModelObjectsTemp
SlabWithToroidalPoreSingletonScore
::do_get_inputs
(Model *m, const ParticleIndexes &pis) const
{
  return IMP::get_particles(m, pis);
}

IMPNPCTRANSPORT_END_NAMESPACE
