/**
 *  \file ExcludeZRangeSingletonScore.h
 *  \brief XXXXXXXXXXXXXX
 *
 *  Copyright 2007-2021 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_EXCLUDE_Z_RANGE_SINGLETON_SCORE_H
#define IMPNPCTRANSPORT_EXCLUDE_Z_RANGE_SINGLETON_SCORE_H

#include "npctransport_config.h"
#include <IMP/SingletonScore.h>
#include <IMP/singleton_macros.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

//! Exclude particles from the given range of z coordinates
class IMPNPCTRANSPORTEXPORT ExcludeZRangeSingletonScore
    : public SingletonScore {
 private:
  double bottom_;  // bottom z coordinated for exclusion
  double top_;     // top z coordinated for exclusion
  double k_;       // the force constant for repulsion out of the z range
 public:
  /**
     Exclude particles from the range of z coordinates [bottom_..top_]
     with repulsive force constant k
   */
  ExcludeZRangeSingletonScore(double bottom, double top, double k);

  /** returns the lowest slab z coordinate */
  double get_bottom_z() const { return bottom_; }

  /** returns the highest slab z coordinate */
  double get_top_z() const { return top_; }

  /** returns the force constant for repulsion out of the z-range */
  double get_k() const { return k_; }

  virtual double evaluate_index(Model *m, ParticleIndex p,
                                DerivativeAccumulator *da) const IMP_OVERRIDE;
  virtual ModelObjectsTemp do_get_inputs(Model *m,
                                         const ParticleIndexes &pis) const
      IMP_OVERRIDE;
  IMP_SINGLETON_SCORE_METHODS(ExcludeZRangeSingletonScore);
  IMP_OBJECT_METHODS(ExcludeZRangeSingletonScore);
};

IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_EXCLUDE_Z_RANGE_SINGLETON_SCORE_H */
