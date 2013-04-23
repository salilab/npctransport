/**
 *  \file ExcludeZRangeSingletonScore.h
 *  \brief XXXXXXXXXXXXXX
 *
 *  Copyright 2007-8 Sali Lab. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_EXCLUDE_Z_RANGE_SINGLETON_SCORE_H
#define IMPNPCTRANSPORT_EXCLUDE_Z_RANGE_SINGLETON_SCORE_H

#include "npctransport_config.h"
#include <IMP/SingletonScore.h>
#include <IMP/singleton_macros.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

class IMPNPCTRANSPORTEXPORT ExcludeZRangeSingletonScore: public SingletonScore
{
 private:
  double bottom_;  // bottom z coordinated for exclusion
  double top_; // top z coordinated for exclusion
  double k_; // the force constant for repulsion out of the z range
 public:
  /**
     Exclude particles from the range of z coordinates [bottom_..top_]
     with repulsive force constant k
   */
  ExcludeZRangeSingletonScore(double bottom, double top, double k);

  IMP_SINGLETON_SCORE(ExcludeZRangeSingletonScore);

  /** returns the lowest slab z coordinate */
  double get_bottom_z() const { return bottom_; }

  /** returns the highest slab z coordinate */
  double get_top_z() const { return top_; }

  /** returns the force constant for repulsion out of the z-range */
  double get_k() const { return k_; }

};

IMPNPCTRANSPORT_END_NAMESPACE

#endif  /* IMPNPCTRANSPORT_EXCLUDE_Z_RANGE_SINGLETON_SCORE_H */
