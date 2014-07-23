/**
 *  \file SitesPairScoreParams.h
 *  \brief A summary of useful information about rigid bodies and their
 *         transformation for eg, caching purposes for SitesPairScore
 *
 *  Copyright 2007-2010 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_INTERNAL_SITES_PAIR_SCORE_PARAMS_H
#define IMPNPCTRANSPORT_INTERNAL_SITES_PAIR_SCORE_PARAMS_H

#include "../npctransport_config.h"

IMPNPCTRANSPORT_BEGIN_INTERNAL_NAMESPACE


/**
   range and k parameters for npctransport::SitesPairScore,
   including some useful cached values
 */
struct SitesPairScoreParams {
  double r; // range
  double k; // interaction k
  double r2; // cached r^2
  double kr; // cachd k*r
  double kr2; // cached k*r^2

  /** initiates structure based on passed range and k params */
  SitesPairScoreParams(double range, double k_coefficient)
  { set_rk(range, k_coefficient); }

  void set_rk(double range, double k_coefficient){
    r = range;
    k = k_coefficient;
    r2 = r*r;
    kr = k*r;
    kr2 = kr*r;
  }

};


IMPNPCTRANSPORT_END_INTERNAL_NAMESPACE

#endif /* IMPNPCTRANSPORT_INTERNAL_SITES_PAIR_SCORE_PARAMS_H */
