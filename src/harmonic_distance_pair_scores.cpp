/**
 *  \file harmonic_distance_pair_scores.cpp
 *  \brief Harmonic scores on the distance between a pair of particles.
 *
 *  Copyright 2007-2019 IMP Inventors. All rights reserved.
 */

#include <IMP/npctransport/harmonic_distance_pair_scores.h>
IMPNPCTRANSPORT_BEGIN_NAMESPACE


HarmonicWellPairScore
::HarmonicWellPairScore
( double rest_length_factor,
  double k, std::string name )
  : PairScore(name),
    rest_length_factor_(rest_length_factor),
    k_(k)
{}

ModelObjectsTemp
HarmonicWellPairScore
::do_get_inputs
( Model *m, const ParticleIndexes &pis) const {
  return IMP::get_particles(m, pis);
}

IMPNPCTRANSPORT_END_NAMESPACE
