/**
 *  \file harmonic_distance_pair_scores.h
 *  \brief Harmonic scores on the distance between a pair of particles.
 *
 *  Copyright 2007-2021 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_HARMONIC_DISTANCE_PAIR_SCORES_H
#define IMPNPCTRANSPORT_HARMONIC_DISTANCE_PAIR_SCORES_H

#include "npctransport_config.h"
#include <IMP/compiler_macros.h>
#include <IMP/PairScore.h>
#include <IMP/pair_macros.h>
#include <IMP/core/XYZR.h>
#include <IMP/algebra/utility.h>

#include <boost/array.hpp>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

/**
   Evaluates a linear pair potential, such that force = k * (delta_length - x0)
   is returned.
   Also, updates the derivative vectors of the particles in the model m.

   @param[out] m    a model for the particles, in which particle derivatives are
                    updated
   @param[in] pp    a pair of particle indices for fast access through internal
                    model methods
   @param[in,out] da accumulator for score derivatives to be updated
   @param[in] delta a vector from pp[1] to pp[0]
   @param delta_length the cached length of delta, assumed correct, and required
                       for faster calculation
   @param x0        resting distance (where score = 0)
   @param k         score linear coefficient. Note that the score is attractive
                    for a positive k, and repulsive for a negative k
                    (assuming the lower the score the better).
   @return the score
 */
inline double do_evaluate_index_harmonic(Model *m, const ParticleIndexPair &pp,
                                DerivativeAccumulator *da,
                                const algebra::Vector3D &delta,
                                double delta_length, double x0, double k) {
  double shifted_length = delta_length - x0;
  double score = 0.5 * k * shifted_length * shifted_length;
  static const double MIN_DISTANCE = .00001;
  if (IMP_LIKELY( da && delta_length > MIN_DISTANCE )) { // Profiling note on use of likely(): in BD simulations, the simulation bottleneck is when da is true, and the spring is likely out of equilibrium
    algebra::Vector3D deriv = (k * shifted_length / delta_length) * delta;  // deriv magnitude is k*shifted_length
    m->add_to_coordinate_derivatives(pp[0], deriv, *da);
    m->add_to_coordinate_derivatives(pp[1], -deriv, *da);
    IMP_LOG(TERSE, "Distance: " << shifted_length << "\nscore: " << score
                                << "\nderiv: " << deriv << std::endl);
  }
  IMP_CHECK_CODE ( else {
      IMP_LOG(TERSE, "Distance: " << shifted_length << "\nscore: " << score
              << std::endl);
    } );
  return score;
}

/**
   A harmonic score between two bonded spheres with a certain
   rest length and spring constant
*/
class IMPNPCTRANSPORTEXPORT HarmonicWellPairScore : public PairScore {
 private:
  double rest_length_factor_, k_;

 public:
  /**
     a linear well pair potential that keeps two particles around
     a resting distance relative to their radii

     @param rest_length_factor the resting distance between particles
                               relative to their sum of radii
     @param k the force constant (attractive if beyond or repulsive
              if below rest length) in units kcal/mol/A^2
     @param name the name of the score
   */
  HarmonicWellPairScore(double rest_length_factor, double k,
                      std::string name = "HarmonicIDPairScore%1%");

  void set_rest_length_factor(double rest_length_factor)
  { rest_length_factor_ = rest_length_factor; }
  double get_rest_length_factor() const { return rest_length_factor_; }
  void set_k(double k)
  { k_ = k; }
  double get_k() { return k_; }
  double evaluate_index(Model *m, const ParticleIndexPair &p,
                        DerivativeAccumulator *da) const IMP_OVERRIDE;
  ModelObjectsTemp do_get_inputs(Model *m,
                  const ParticleIndexes &pis) const IMP_OVERRIDE;
  IMP_PAIR_SCORE_METHODS(HarmonicWellPairScore);
  IMP_OBJECT_METHODS(HarmonicWellPairScore);
  ;
};

#ifndef IMP_DOXYGEN
inline double HarmonicWellPairScore::evaluate_index(
    Model *m, const ParticleIndexPair &pp, DerivativeAccumulator *da) const {
  IMP_OBJECT_LOG;

  algebra::Sphere3D s0 = m->get_sphere(pp[0]);
  algebra::Sphere3D s1 = m->get_sphere(pp[1]);
  double x0 = (s0.get_radius() + s1.get_radius()) * rest_length_factor_;
  algebra::Vector3D delta = s0.get_center() - s1.get_center();
  double delta_length_2 = delta.get_squared_magnitude();
  double delta_length = std::sqrt(delta_length_2);
  return  // k_ > 0 = get spheres closer
    do_evaluate_index_harmonic(m, pp, da, delta, delta_length, x0, k_);
}
#endif

IMP_OBJECTS(HarmonicWellPairScore, HarmonicWellPairScores);


IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_HARMONIC_DISTANCE_PAIR_SCORES_H */
