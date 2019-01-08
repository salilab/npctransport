/**
 *  \file linear_distance_pair_scores.cpp
 *  \brief Linear scores on the distance between a pair of particles.
 *
 *  Copyright 2007-2019 IMP Inventors. All rights reserved.
 */

#include <IMP/npctransport/linear_distance_pair_scores.h>
IMPNPCTRANSPORT_BEGIN_NAMESPACE

LinearSoftSpherePairScore::LinearSoftSpherePairScore(double k, std::string name)
    : PairScore(name), k_(k) {}

ModelObjectsTemp LinearSoftSpherePairScore::do_get_inputs(
    Model *m, const ParticleIndexes &pis) const {
  return IMP::get_particles(m, pis);
}

LinearInteractionPairScore::LinearInteractionPairScore(double k_rep,
                                                       double range_attr,
                                                       double k_attr,
                                                       std::string name)
    : PairScore(name), range_attr_(range_attr), k_rep_(k_rep), k_attr_(k_attr) {}

ModelObjectsTemp LinearInteractionPairScore::do_get_inputs(
    Model *m, const ParticleIndexes &pis) const {
  return IMP::get_particles(m, pis);
}

LinearWellPairScore::LinearWellPairScore(double rest_length_factor,
                                         double k, std::string name)
    : PairScore(name), rest_length_factor_(rest_length_factor), k_(k)
{}

ModelObjectsTemp LinearWellPairScore::do_get_inputs(
    Model *m, const ParticleIndexes &pis) const {
  return IMP::get_particles(m, pis);
}

IMPNPCTRANSPORT_END_NAMESPACE
