/**
 *  \file functor_linear_distance_pair_scores.h
 *  \brief A Score on the distance between a pair of particles.
 *
 *  Copyright 2007-2022 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_FUNCTOR_LINEAR_DISTANCE_PAIR_SCORES_H
#define IMPNPCTRANSPORT_FUNCTOR_LINEAR_DISTANCE_PAIR_SCORES_H

#include "npctransport_config.h"
#include <IMP/npctransport/functor_linear_distance_pair_scores_typedefs.h>
#include <IMP/algebra/utility.h>
#include <IMP/score_functor/LinearLowerBound.h>
#include <IMP/score_functor/SphereDistance.h>
#include <IMP/score_functor/DistancePairScore.h>
#include <IMP/score_functor/distance_pair_score_macros.h>

#include <boost/array.hpp>

IMPNPCTRANSPORT_BEGIN_NAMESPACE


//! Score a pair of particles
class FunctorLinearSoftSpherePairScore
    : public
//#ifndef SWIG
      IMP::score_functor::DistancePairScore<LinearSoftSphereScore>
//#else
//     IMP::PairScore
//#endif
      {
  typedef IMP::score_functor::DistancePairScore<LinearSoftSphereScore> P;

 public:
  FunctorLinearSoftSpherePairScore(double k,
                                   std::string name = "LinearSSPairScore%1%")
      : P(LinearSoftSphereScore(k), name) {}
  //#ifdef SWIG
  IMP_OBJECT_METHODS(FunctorLinearSoftSpherePairScore);
  //#endif
};


typedef score_functor::SphereDistance<LinearInteraction> LinearInteractionScore;

//! Score a pair of particles
class FunctorLinearInteractionPairScore
    : public
//#ifndef SWIG
      IMP::score_functor::DistancePairScore<LinearInteractionScore>
//#else
//      IMP::PairScore
//#endif
      {
  typedef IMP::score_functor::DistancePairScore<LinearInteractionScore> P;

 public:
  FunctorLinearInteractionPairScore(double krep, double attr_range,
                                    double kattr,
                                    std::string name = "LinearSSPairScore%1%")
      : P(LinearInteractionScore(LinearInteraction(krep, attr_range, kattr)),
          name) {}
  //#ifdef SWIG
  IMP_OBJECT_METHODS(FunctorLinearInteractionPairScore);
  //#endif
};

IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_FUNCTOR_LINEAR_DISTANCE_PAIR_SCORES_H */
