/**
 *  \file HarmonicSpringSingletonScore.cpp
 *  \brief Harmonic scores on two particles tethered to a spring
 *
 *  Copyright 2007-2021 IMP Inventors. All rights reserved.
 */

#include <IMP/npctransport/HarmonicSpringSingletonScore.h>
#include <IMP/npctransport/RelaxingSpring.h>
IMPNPCTRANSPORT_BEGIN_NAMESPACE


HarmonicSpringSingletonScore
::HarmonicSpringSingletonScore
( double k1,
  double k2,
  std::string name )
  : SingletonScore(name),
    k1_(k1),
    k2_(k2)
{}

ModelObjectsTemp
HarmonicSpringSingletonScore
::do_get_inputs
( Model *m, const ParticleIndexes &pis) const {
  ModelObjectsTemp ret(3 * pis.size());
  for (unsigned int i = 0; i < pis.size(); ++i) {
    RelaxingSpring s(m, pis[i]);
    ret[3 * i + 0] = s.get_bonded_particle_0();
    ret[3 * i + 1] = s.get_bonded_particle_1();
    ret[3 * i + 2] = s.get_particle();
  }
  return ret;
}

IMPNPCTRANSPORT_END_NAMESPACE
