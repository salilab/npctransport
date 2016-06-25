/**
 *  \file SitesPairScoreParameters.h
 *  \brief A summary of useful information about rigid bodies and their
 *         transformation for eg, caching purposes for SitesPairScore
 *
 *  Copyright 2007-2010 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_SITES_PAIR_SCORE_PARAMS_H
#define IMPNPCTRANSPORT_SITES_PAIR_SCORE_PARAMS_H

#include "npctransport_config.h"
#include <IMP/value_macros.h>
#include <IMP/showable_macros.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE


/**
   range and k parameters for npctransport::SitesPairScore,
   including some useful cached values
 */
struct IMPNPCTRANSPORTEXPORT SitesPairScoreParameters {
  double r; // range
  double k; // interaction k
  double r2; // cached r^2
  double kr; // cachd k*r
  double kr2; // cached k*r^2
  double cosSigma1_max;
  double cosSigma2_max;
  bool is_orientational;

  /** initiates structure based on passed range and k params 

      @param range range of attraction in A
      @param k_coefficient k in units of [kCal/mol/A] or [kCal/mol/A^2] (for isotropic and anisotropic, respectively - isotropic when either sigma_max is 0.0)
      @param sigma1_max maximal angle for site 1 [0.0 if isotropic interaction]
      @param sigma2_max maximal angle for site 2 [0.0 if isotropic interaction]
   */
  SitesPairScoreParameters(double range, double k_coefficient, double sigma1_max_deg=0.0, double sigma2_max_deg=0.0)
  { 
    set_range(range);
    set_force_coefficient(k_coefficient);
    set_sigma1_max(sigma1_max_deg);
    set_sigma2_max(sigma2_max_deg);
    update_is_orientational();
  }

  void set_range(double r){
    this->r=r;
    update_cache();
  }

  void set_force_coefficient(double k){
    this->k=k;
    update_cache();
  }

  void set_sigma1_max(double sigma1_max_deg){
    double sigma1_max_rad= sigma1_max_deg * 4.0 * atan (1.0) / 180.0;
    cosSigma1_max=std::cos(sigma1_max_rad);
    update_is_orientational();
  }

  void set_sigma2_max(double sigma2_max_deg){
    double sigma2_max_rad= sigma2_max_deg * 4.0 * atan (1.0) / 180.0;
    cosSigma2_max=std::cos(sigma2_max_rad);
    update_is_orientational();
  }

  void update_is_orientational(){
    is_orientational= std::abs((cosSigma1_max-1.0) + (cosSigma2_max-1.0)) < 0.0001;
  }

  private:
    void update_cache(){
      r2 = r*r;
      kr = k*r;
      kr2 = kr*r;
    }

 public:
    IMP_SHOWABLE_INLINE(SitesPairScoreParameters, 
		      out << "sites pair score params"
		      << " range " << r 
		      << " k " << k
		      << " cos(sigma1_max) " << cosSigma1_max   
		      << " cos(sigma2_max) " << cosSigma2_max   
		      << " is_orientational " << is_orientational
		      << std::endl);


};

IMP_VALUES(SitesPairScoreParameters,SitesPairScoreParametersList);

IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_SITES_PAIR_SCORE_PARAMS_H */
