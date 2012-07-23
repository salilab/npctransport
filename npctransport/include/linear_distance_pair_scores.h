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

   @param m[out]    a model for the particles, in which particle derivatives are
                    updated
   @param pp[in]    a pair of particle indices for fast access through internal
                    model methods
   @param delta[in] a vector that represents the displacement between the
                    two particles
   @param x0        resting distance (where score = 0)
   @param k         score linear coefficient. Note that the score is attractive
                    for a positive k, and repulsive for a negative k
                    (assuming the lower the score the better).
   @return the score
 */
inline double do_evaluate_index(Model *m, const ParticleIndexPair& pp,
                                DerivativeAccumulator *da,
                                const algebra::Vector3D &delta,
                                double delta_length,
                                double x0, double k) {
  double shifted_length=delta_length-x0;
  double score= k*shifted_length;
  static const double MIN_DISTANCE = .00001;
  if (da && delta_length > MIN_DISTANCE) {
    algebra::Vector3D deriv= k*delta/delta_length; // deriv magnitude is k
    m->add_to_coordinate_derivatives(pp[0], deriv, *da);
    m->add_to_coordinate_derivatives(pp[1], -deriv, *da);
    IMP_LOG(TERSE, "Distance: " << shifted_length
            << "\nscore: " << score << "\nderiv: " << deriv << std::endl);
  } else {
    IMP_LOG(TERSE, "Distance: " << shifted_length
            << "\nscore: " << score << std::endl);
  }
  return score;
}

/**
   A soft linear repulsive score between two spheres.
   If the spheres penetrate each other, the score is
   ( k_ * [penetration-magnitude] ), otherwise the score is 0.
 */
class IMPNPCTRANSPORTEXPORT LinearSoftSpherePairScore:
    public PairScore {
    private:
      // score potential coefficient (larger = stronger repulsion)
      double k_;

    public:
      LinearSoftSpherePairScore(double k,
                                std::string name="LinearSSPairScore%1%"
                                );
      IMP_INDEX_PAIR_SCORE(LinearSoftSpherePairScore);
    };


// evaluates the soft linear repulsive score, using the internal Model
// indexing to access particles data (see Model/include/internal/)
inline double LinearSoftSpherePairScore
::evaluate_index(Model *m, const ParticleIndexPair& pp,
                 DerivativeAccumulator *da) const {
  IMP_OBJECT_LOG;
  // compuate distances into cached variables
  algebra::Vector3D delta=
    m->get_sphere(pp[0]).get_center()
    - m->get_sphere(pp[1]).get_center();
  double delta_length_2=
    delta.get_squared_magnitude();   // (work with sqr for efficiency)
  double x0= m->get_sphere(pp[0]).get_radius()
    + m->get_sphere(pp[1]).get_radius();
  double x0_2 = x0 * x0;
  bool not_penetrating = delta_length_2 > x0_2;
  if ( not_penetrating )
    return 0;
  // penetrating spheres ; -k_ in order to get a repulsive score
  double delta_length = std::sqrt(delta_length_2);
  return do_evaluate_index
    (m, pp, da, delta, delta_length, x0, -k_);
}




/**
   A soft linear attractive / repulsive score between two spheres.
   The score is 0 if the spheres are beyond the attractive range.
   Within the attractive range, the score decreases linearly (= attraction)
   with slope k_attr_ until the spheres touch. Once the spheres begin to
   penetrate each other, the score rises linearly with slope k_rep_
   (= repulsion), though it may be negative for small penetration.
*/
class IMPNPCTRANSPORTEXPORT LinearInteractionPairScore:
    public PairScore
{
#ifndef SWIG
 public:
  /** cache for reusable intermediate calculation */
  struct EvaluationCache{
    // the squared-distance between particle centers:
    double particles_delta_squared;
    // the sum of particle radii:
    double sum_particles_radii;
  };
#endif

 private:
  double
    attr_range_, // range of attraction between particles
    k_rep_,      // repulsion coeficcient
    k_attr_;     // attraction coefficient

  // cache for reusable intermediate calcualtions of evaluate_infex
  mutable EvaluationCache cache_;

public:
  LinearInteractionPairScore(double krep, double attr_range, double kattr,
                             std::string name="LinearIDPairScore%1%");

#ifndef SWIG
  // returns cached intermediate computations from last call to
  // evaluate_index()
  EvaluationCache const& get_evaluation_cache() const
  { return cache_; }
#endif


  IMP_INDEX_PAIR_SCORE(LinearInteractionPairScore);
};

inline double LinearInteractionPairScore
::evaluate_index(Model *m, const ParticleIndexPair& pp,
                 DerivativeAccumulator *da) const {
  IMP_OBJECT_LOG;

  // Associate intermediate variables with cache_, for further reuse:
  double& delta_length_2 = cache_.particles_delta_squared;
  double& x0 =             cache_.sum_particles_radii;
  algebra::Vector3D delta =
    m->get_sphere(pp[0]).get_center()
    - m->get_sphere(pp[1]).get_center();
  delta_length_2 =
    delta.get_squared_magnitude();
  x0 =
    m->get_sphere(pp[0]).get_radius()
    + m->get_sphere(pp[1]).get_radius();
  // Terminate immediately if very far, work with squares for speed
  // equivalent to [delta_length > x0 + attr_range]:
  if ( delta_length_2 > std::pow(x0 + attr_range_, 2) )
    return 0;
  double offset = -attr_range_*k_attr_;
  double delta_length = std::sqrt(delta_length_2);
  if (delta_length > x0) { // attractive regime
    return // decreases with slope k_attr_ as spheres get closer
      ( do_evaluate_index(m, pp, da, delta, delta_length, x0, k_attr_)
        + offset );
  } else { // repulsive regime, may be negative for small penetration
    return
      ( do_evaluate_index(m, pp, da, delta, delta_length, x0, -k_rep_)
        + offset );
  }
}






/**
   A soft linear attractive / repulsive score between two spheres
   for the backbone of a chain.
*/
class IMPNPCTRANSPORTEXPORT LinearWellPairScore:
    public PairScore
{
 private:
  double x0_, k_;
public:
  LinearWellPairScore(double x0, double k,
                             std::string name="LinearIDPairScore%1%");

  void set_x0(double x0) {
    x0_=x0;
  }
  double get_x0() const {
    return x0_;
  }
  IMP_INDEX_PAIR_SCORE(LinearWellPairScore);
};

inline double LinearWellPairScore
::evaluate_index(Model *m, const ParticleIndexPair& pp,
                 DerivativeAccumulator *da) const {
  IMP_OBJECT_LOG;

  // Associate intermediate variables with cache_, for further reuse:
  algebra::Vector3D delta =  m->get_sphere(pp[0]).get_center()
    - m->get_sphere(pp[1]).get_center();
  double delta_length_2 =
    delta.get_squared_magnitude();
  double delta_length= (1.0)/(1.0/std::sqrt(delta_length_2));
  if (delta_length > x0_) { // attractive regime
    return // decreases with slope k_attr_ as spheres get closer
        do_evaluate_index(m, pp, da, delta, delta_length, x0_, k_);
  } else { // repulsive regime, may be negative for small penetration
    return
      do_evaluate_index(m, pp, da, delta, delta_length, x0_, -k_);
  }
}

IMP_OBJECTS(LinearWellPairScore, LinearWellPairScores);

IMPNPCTRANSPORT_END_NAMESPACE

#endif  /* IMPNPCTRANSPORT_LINEAR_DISTANCE_PAIR_SCORES_H */
