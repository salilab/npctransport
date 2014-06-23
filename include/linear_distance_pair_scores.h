/**
 *  \file linear_distance_pair_scores.h
 *  \brief A Score on the distance between a pair of particles.
 *
 *  Copyright 2007-2010 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_LINEAR_DISTANCE_PAIR_SCORES_H
#define IMPNPCTRANSPORT_LINEAR_DISTANCE_PAIR_SCORES_H

#include "npctransport_config.h"
#include <IMP/PairScore.h>
#include <IMP/pair_macros.h>
#include <IMP/core/XYZR.h>
#include <IMP/algebra/utility.h>
#include "internal/sites.h"

#include <boost/array.hpp>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

/**
   Evaluates a linear pair potential, such that score = k * (delta_length - x0)
   is returned.
   Also, updates the derivative vectors of the particles in the model m.

   @param[out] m    a model for the particles, in which particle derivatives are
                    updated
   @param[in] pp    a pair of particle indices for fast access through internal
                    model methods
   @param[in,out] da accumulator for score derivatives to be updated
   @param[in] delta a vector that represents the displacement between the
                    two particles
   @param delta_length the cached length of delta, assumed correct, and required
                       for faster calculation
   @param x0        resting distance (where score = 0)
   @param k         score linear coefficient. Note that the score is attractive
                    for a positive k, and repulsive for a negative k
                    (assuming the lower the score the better).
   @return the score
 */
inline double do_evaluate_index(Model *m, const ParticleIndexPair &pp,
                                DerivativeAccumulator *da,
                                const algebra::Vector3D &delta,
                                double delta_length, double x0, double k) {
  double shifted_length = delta_length - x0;
  double score = k * shifted_length;
  static const double MIN_DISTANCE = .00001;
  if (da && delta_length > MIN_DISTANCE) {
    algebra::Vector3D deriv = k * delta / delta_length;  // deriv magnitude is k
    m->add_to_coordinate_derivatives(pp[0], deriv, *da);
    m->add_to_coordinate_derivatives(pp[1], -deriv, *da);
    IMP_LOG(TERSE, "Distance: " << shifted_length << "\nscore: " << score
                                << "\nderiv: " << deriv << std::endl);
  } else {
    IMP_LOG(TERSE, "Distance: " << shifted_length << "\nscore: " << score
                                << std::endl);
  }
  return score;
}

/**
   A soft linear repulsive score between two spheres.
   If the spheres penetrate each other, the score is
   ( k_ * [penetration-magnitude] ), otherwise the score is 0.
 */
class IMPNPCTRANSPORTEXPORT LinearSoftSpherePairScore : public PairScore {
 private:
  // score potential coefficient (larger = stronger repulsion)
  double k_;

 public:

  LinearSoftSpherePairScore(double k,
                            std::string name = "LinearSSPairScore%1%");

  virtual double evaluate_index(Model *m, const ParticleIndexPair &p,
                                DerivativeAccumulator *da) const IMP_OVERRIDE;

  virtual ModelObjectsTemp do_get_inputs(Model *m,
                                         const ParticleIndexes &pis) const;

  //! returns the k for sphere-sphere repulsion
  double get_k() const { return k_; }

  IMP_PAIR_SCORE_METHODS(LinearSoftSpherePairScore);
  IMP_OBJECT_METHODS(LinearSoftSpherePairScore);
  ;
};

#ifndef IMP_DOXYGEN
// evaluates the soft linear repulsive score, using the internal Model
// indexing to access particles data (see Model/include/internal/)
inline double LinearSoftSpherePairScore::evaluate_index(
    Model *m, const ParticleIndexPair &pp, DerivativeAccumulator *da) const {
  IMP_OBJECT_LOG;
  // compuate distances into cached variables
  algebra::Vector3D delta =
      m->get_sphere(pp[0]).get_center() - m->get_sphere(pp[1]).get_center();
  double delta_length_2 =
      delta.get_squared_magnitude();  // (work with sqr for efficiency)
  double x0 =
      m->get_sphere(pp[0]).get_radius() + m->get_sphere(pp[1]).get_radius();
  double x0_2 = x0 * x0;
  bool not_penetrating = delta_length_2 > x0_2;
  if (not_penetrating) return 0;
  // penetrating spheres ; -k_ in order to get a repulsive score
  double delta_length = std::sqrt(delta_length_2);
  return do_evaluate_index(m, pp, da, delta, delta_length, x0, -k_);
}
#endif

/**
   A soft linear attractive / repulsive score between two spheres.
   The score is 0 if the spheres are beyond the attractive range.
   Within the attractive range, the score decreases linearly (= attraction)
   with slope k_attr_ until the spheres touch. Once the spheres begin to
   penetrate each other, the score rises linearly with slope k_rep_
   (= repulsion), though it may be negative for small penetration.
*/
class IMPNPCTRANSPORTEXPORT LinearInteractionPairScore : public PairScore {
#ifndef SWIG
 public:
  /** cache for reusable intermediate calculation */
  struct EvaluationCache {
    // the squared-distance between particle centers:
    double particles_delta_squared;
    // the sum of particle radii:
    double sum_particles_radii;
  };
#endif

 private:
  double range_attr_,  // range of attraction between particles
      k_rep_,          // repulsion coeficcient
      k_attr_;         // attraction coefficient

  // cache for reusable intermediate calcualtions of evaluate_infex
  mutable EvaluationCache cache_;

 public:
  LinearInteractionPairScore(double k_rep, double range_attr, double k_attr,
                             std::string name = "LinearIDPairScore%1%");

#ifndef SWIG
  // returns cached intermediate computations from last call to
  // this->evaluate_index()
  EvaluationCache const &get_evaluation_cache() const { return cache_; }
#endif

  IMP_IMPLEMENT_INLINE(double evaluate(const ParticlePair &p,
                                       DerivativeAccumulator *da) const,
                       {
    return evaluate_index(IMP::internal::get_model(p),
                          IMP::internal::get_index(p), da);
  });
  IMP_IMPLEMENT(double evaluate_index(Model *m, const ParticleIndexPair &p,
                                      DerivativeAccumulator *da) const
                IMP_OVERRIDE);
  IMP_IMPLEMENT_INLINE(double evaluate_if_good_index(Model *m,
                                                     const ParticleIndexPair &p,
                                                     DerivativeAccumulator *da,
                                                     double max) const,
                       {
    IMP_UNUSED(max);
    return evaluate_index(m, p, da);
  });
  virtual double evaluate_indexes(Model *m, const ParticleIndexPairs &p,
                                  DerivativeAccumulator *da,
                                  unsigned int lower_bound,
                                  unsigned int upper_bound) const IMP_OVERRIDE {
    double ret = 0;
    for (unsigned int i = lower_bound; i < upper_bound; ++i) {
      ret += evaluate_index(m, p[i], da);
    }
    return ret;
  }
  double evaluate_if_good_index(Model *m, const ParticleIndexPairs &p,
                                DerivativeAccumulator *da, double max,
                                unsigned int lower_bound,
                                unsigned int upper_bound) const {
    double ret = 0;
    for (unsigned int i = lower_bound; i < upper_bound; ++i) {
      ret += evaluate_if_good_index(m, p[i], da, max - ret);
      if (ret > max) return std::numeric_limits<double>::max();
    }
    return ret;
  }
  ModelObjectsTemp do_get_inputs(Model *m, const ParticleIndexes &pis) const
      IMP_OVERRIDE;

  //! returns the range for sphere-sphere attraction
  double get_range_attraction() const { return range_attr_; }

  //! returns the k for sphere-sphere attraction
  double get_k_attraction() const { return k_attr_; }

  //! returns the k for sphere-sphere repulsion
  double get_k_repulsion() const { return k_rep_; }


  IMP_OBJECT_METHODS(LinearInteractionPairScore);
};

#ifndef IMP_DOXYGEN
inline double LinearInteractionPairScore::evaluate_index(
    Model *m, const ParticleIndexPair &pp, DerivativeAccumulator *da) const {
  IMP_OBJECT_LOG;

  // Associate intermediate variables with cache_, for further reuse:
  double &delta_length_2 = cache_.particles_delta_squared;
  double &x0 = cache_.sum_particles_radii;
  algebra::Vector3D delta =
      m->get_sphere(pp[0]).get_center() - m->get_sphere(pp[1]).get_center();
  delta_length_2 = delta.get_squared_magnitude();
  x0 = m->get_sphere(pp[0]).get_radius() + m->get_sphere(pp[1]).get_radius();
  // Terminate immediately if very far, work with squares for speed
  // equivalent to [delta_length > x0 + attr_range]:
  if (delta_length_2 > std::pow(x0 + range_attr_, 2)) return 0;
  double offset = -range_attr_ * k_attr_;
  double delta_length = std::sqrt(delta_length_2);
  if (delta_length > x0) {  // attractive regime
    return  // decreases with slope k_attr_ as spheres get closer
        (do_evaluate_index(m, pp, da, delta, delta_length, x0, k_attr_) +
         offset);
  } else {  // repulsive regime, may be negative for small penetration
    return (do_evaluate_index(m, pp, da, delta, delta_length, x0, -k_rep_) +
            offset);
  }
}

#endif

/**
   A soft linear attractive / repulsive score between two spheres
   for the backbone of a chain.
*/
class IMPNPCTRANSPORTEXPORT LinearWellPairScore : public PairScore {
 private:
  double rest_length_factor_, k_;

 public:
  /**
     a linear well pair potential that keeps two particles around
     a resting distance relative to their radii

     @param rest_length_factor the resting distance between particles
                               relative to their sum of radii
     @param k the force constant (attractive if beyond or repulsive
              if below rest length)
     @param name the name of the score
   */
  LinearWellPairScore(double rest_length_factor, double k,
                      std::string name = "LinearIDPairScore%1%");

  void set_rest_length_factor(double rest_length_factor)
  { rest_length_factor_ = rest_length_factor; }
  double get_rest_length_factor() const { return rest_length_factor_; }
  void set_k(double k)
  { k_ = k; }
  double get_k() { return k_; }
  double evaluate_index(Model *m, const ParticleIndexPair &p,
                        DerivativeAccumulator *da) const IMP_OVERRIDE;
  ModelObjectsTemp do_get_inputs(Model *m, const ParticleIndexes &pis) const;
  IMP_PAIR_SCORE_METHODS(LinearWellPairScore);
  IMP_OBJECT_METHODS(LinearWellPairScore);
  ;
};

#ifndef IMP_DOXYGEN
inline double LinearWellPairScore::evaluate_index(
    Model *m, const ParticleIndexPair &pp, DerivativeAccumulator *da) const {
  IMP_OBJECT_LOG;

  algebra::Sphere3D s0 = m->get_sphere(pp[0]);
  algebra::Sphere3D s1 = m->get_sphere(pp[1]);
  double x0 = (s0.get_radius() + s1.get_radius()) * rest_length_factor_;
  algebra::Vector3D delta = s0.get_center() - s1.get_center();
  double delta_length_2 = delta.get_squared_magnitude();
  double delta_length = std::sqrt(delta_length_2);
  if (delta_length > x0) {  // attractive regime
    return  // k_ > 0 = get spheres closer
        do_evaluate_index(m, pp, da, delta, delta_length, x0, k_);
  } else {  // -k_ < 0 = keep spheres apart
    return do_evaluate_index(m, pp, da, delta, delta_length, x0, -k_);
  }
}
#endif

IMP_OBJECTS(LinearWellPairScore, LinearWellPairScores);


IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_LINEAR_DISTANCE_PAIR_SCORES_H */
