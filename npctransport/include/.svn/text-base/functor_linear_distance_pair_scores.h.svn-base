/**
 *  \file functor_linear_distance_pair_scores.h
 *  \brief A Score on the distance between a pair of particles.
 *
 *  Copyright 2007-2010 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_FUNCTOR_LINEAR_DISTANCE_PAIR_SCORES_H
#define IMPNPCTRANSPORT_FUNCTOR_LINEAR_DISTANCE_PAIR_SCORES_H

#include "npctransport_config.h"
#include <IMP/algebra/utility.h>
#include <IMP/score_functor/LinearLowerBound.h>
#include <IMP/score_functor/SphereDistance.h>
#include <IMP/score_functor/DistancePairScore.h>
#include <IMP/score_functor/distance_pair_score_macros.h>

#include <boost/array.hpp>



IMPNPCTRANSPORT_BEGIN_NAMESPACE


typedef score_functor::SphereDistance<score_functor::LinearLowerBound>
LinearSoftSphereScore;
IMP_FUNCTOR_DISTANCE_PAIR_SCORE(FunctorLinearSoftSpherePairScore,
                                LinearSoftSphereScore,
                                (double k,
                                 std::string name="LinearSSPairScore%1%"), (k));


/**
   A soft linear attractive / repulsive score between two spheres.
   The score is 0 if the spheres are beyond the attractive range.
   Within the attractive range, the score decreases linearly (= attraction)
   with slope k_attr_ until the spheres touch. Once the spheres begin to
   penetrate each other, the score rises linearly with slope k_rep_
   (= repulsion), though it may be negative for small penetration.
*/
class LinearInteraction: public score_functor::LinearLowerBound {
  typedef score_functor::LinearLowerBound P;
  double attr_range_; // range of attraction between particles
  double k_attr_;     // attraction coefficient
public:
  LinearInteraction(double krep, double attr_range, double kattr):
    P(krep) {
    attr_range_=attr_range;
    k_attr_=kattr;
  }
  // depend on get_is_trivially_zero
  template <unsigned int D>
  double get_score(Model *m, const ParticleIndexTuple<D>&pp,
                   double distance) const {
    if (distance < 0) {
      return P::get_score(m, pp, distance)-k_attr_*attr_range_;
    } else {
      IMP_USAGE_CHECK(distance <= attr_range_,
                      "It is trivially 0.");
      return k_attr_*(distance-attr_range_);
    }
  }
  template <unsigned int D>
  DerivativePair get_score_and_derivative(Model *m,
                                          const ParticleIndexTuple<D>&p,
                                          double distance) const {
    if (distance< 0) {
      DerivativePair dp= P::get_score_and_derivative(m, p, distance);
      return DerivativePair(dp.first-k_attr_*attr_range_,
                            dp.second);
    } else  {
      return DerivativePair(k_attr_*(distance-attr_range_),
                            k_attr_);
    }
  }
  template <unsigned int D>
  double get_maximum_range(Model *, const ParticleIndexTuple<D>& ) const {
    return attr_range_;
  }
  template <unsigned int D>
  bool get_is_trivially_zero(Model *, const ParticleIndexTuple<D>& ,
                             double squared_distance) const {
    return squared_distance > algebra::get_squared(attr_range_);
  }
};


typedef score_functor::SphereDistance<LinearInteraction>
LinearInteractionScore;
IMP_FUNCTOR_DISTANCE_PAIR_SCORE(FunctorLinearInteractionPairScore,
                                LinearInteractionScore,
                                (double krep, double attr_range, double kattr,
                                 std::string name="LinearSSPairScore%1%"),
                                (LinearInteraction(krep, attr_range, kattr)));


IMPNPCTRANSPORT_END_NAMESPACE

#endif  /* IMPNPCTRANSPORT_FUNCTOR_LINEAR_DISTANCE_PAIR_SCORES_H */
