/**
 *  \file DistancePairScore.cpp
 *  \brief A Score on the distance between a pair of particles.
 *
 *  Copyright 2007-2010 IMP Inventors. All rights reserved.
 */

#include <IMP/npctransport/SitesPairScore.h>
#include <IMP/core/XYZR.h>
#include <IMP/core/rigid_bodies.h>
#include <IMP/generic.h>
#include <boost/preprocessor/seq/for_each_product.hpp>
#include <boost/preprocessor/seq/to_tuple.hpp>
#include <boost/preprocessor/facilities/apply.hpp>
#include <boost/preprocessor/seq/elem.hpp>
#include <boost/preprocessor/tuple/to_seq.hpp>
#include <math.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

namespace {
  // return maximum sum of distance between pairs of particles with
  // site list R0 and R1 respectively, such that there exists a site
  // r0 in R0 and a site r1 in R1 that overlap. There are the pair of
  // sites that has the maximal sum of site radii and distances from
  // the center.
  double get_max_r_sum(const algebra::Sphere3Ds& R0,
                const algebra::Sphere3Ds& R1) {
    double max = 0.0;
    for(unsigned int i = 0; i < R0.size(); i++){
      for(unsigned int j = 0; j < R1.size(); j++){
        double s = (R0[i].get_center().get_magnitude() +  R0[i].get_radius())
          + (R1[j].get_center().get_magnitude() +  R1[j].get_radius());
        max = std::max(max,s);
      }
    }
    return max;
  }
}


// orientation-dependent interaction score
SitesPairScore::SitesPairScore(double range, double k,
			       double sigma0_deg, double sigma1_deg,
                               double range_nonspec_attraction,
                               double k_nonspec_attraction,
                               double k_repulsion,
                               const algebra::Sphere3Ds &sites0,
                               const algebra::Sphere3Ds &sites1)
  :
  P(k_repulsion, range_nonspec_attraction, k_nonspec_attraction,
    "SitesPairScore %1%"),
  params_(range, k, sigma0_deg, sigma1_deg),
  sites0_(sites0),
  sites1_(sites1),
  is_cache_active_(false),
  cur_cache_id_(INVALID_CACHE_ID)
{
  IMP_LOG_PROGRESS( "Setting up SitesPairScore with sites0 "
		    << sites0_ << " sites1 " << sites1_ << std::endl);
  is_orientational_score_ = (sigma0_deg > 0.0 && sigma1_deg > 0.0);
  if(!is_orientational_score_){
    IMP_LOG(WARNING, "Creating old version of SitesPairScore, for "
	    "backwards compatibility" << std::endl);
  }
  // Find upper bound for distance between particles whose sites interact
  // to be used for fast filtering - the range + sites radii
  double ubound_distance = (get_max_r_sum(sites0, sites1) + params_.r);
  ubound_distance2_ = ubound_distance * ubound_distance;
  IMP_LOG_PROGRESS( "UBOUND,2: " << ubound_distance << " ; "
                    << this->ubound_distance2_ << std::endl);
}

ModelObjectsTemp SitesPairScore::do_get_inputs(
    Model *m, const ParticleIndexes &pis) const {
  return IMP::get_particles(m, pis);
}

// the sites of each particle are transformed to a common frame of reference
// (using the reference frame of each particle), and the site-specific
// attraction, and the inter-particle non specific attraction and repulsion
// are evaluated and summed.
inline double SitesPairScore::evaluate_index(Model *m,
                                             const ParticleIndexPair &pip,
                                             DerivativeAccumulator *da) const {
  IMP_OBJECT_LOG;

  // I. evaluate non-specific attraction and repulsion between
  //    parent particles before computing for specific sites :
  double non_specific_score = P::evaluate_index(m, pip, da);
  LinearInteractionPairScore::EvaluationCache const&
    lips_cache= P::get_evaluation_cache();

  // II. Return if parent particles are out of site-specific interaction range
  //     using cache to avoid some redundant calcs
  double const& distance2= lips_cache.particles_delta_squared;
  IMP_LOG(PROGRESS, "distance2 " << distance2
          << " ; distance upper-bound " << ubound_distance2_ <<  std::endl);
  if (distance2 > ubound_distance2_) {
    IMP_LOG(PROGRESS, "Sites contribution is 0.0 and non-specific score is "
            << non_specific_score << std::endl);
    return non_specific_score;
  }

  // III. evaluate site-specific contributions :
  // bring sites_ to the frame of reference of nn_sites_ and nn_
  ParticleIndex pi0 = pip[0];
  ParticleIndex pi1 = pip[1];
  // get rbi0/1 from cache, update if needed
  internal::RigidBodyInfo rbi0 = get_rigid_body_info(m, pi0);
  internal::RigidBodyInfo rbi1 = get_rigid_body_info(m, pi1);
  IMP_LOG_PROGRESS( "RBI0.rb " << rbi0.rb
                    << " RB0.cache_id " << rbi0.cache_id
                    << "RBI0.tr " << rbi0.tr << std::endl);
  IMP_LOG_PROGRESS( "RBI1.cache_id " << rbi1.cache_id
                    << "RBI1.RB " << rbi1.rb
                    << "RBI1.tr " << rbi1.tr << std::endl);
  // sum over specific interactions between all pairs of sites:
  double sum = 0;
  for (unsigned int i = 0; i < sites0_.size(); ++i) {
    algebra::Vector3D g0 = rbi0.tr.get_transformed(sites0_[i].get_center());
    for(unsigned int j = 0 ; j < sites1_.size(); ++j) {
      algebra::Vector3D g1 = rbi1.tr.get_transformed(sites1_[j].get_center());
      IMP_LOG_PROGRESS( "Evaluating sites at global coordinates: " << g0 << " ; " << g1 << std::endl );
      if(is_orientational_score_)
        {
          sum +=
            internal::evaluate_pair_of_sites(params_,
                                          rbi0, rbi1,
                                          sites0_[i], sites1_[j],
                                          g0, g1,
                                          da);
      } else
        { // old score
          sum +=
            internal::evaluate_one_site_3(params_.k,
                                          params_.r,
                                          rbi0, rbi1,
                                          sites0_[i], sites1_[j],
                                          g0, g1,
                                          da);
        }
      IMP_LOG_PROGRESS( "Sum " << sum << std::endl);
    }
  }
  IMP_LOG(PROGRESS,
          "Sites contribution is " << sum <<
          " and non-specific soft sphere contribution is " << non_specific_score
          << std::endl);
  return sum + non_specific_score;
}

// Restraints SitesPairScore::do_create_current_decomposition(
//     Model *m, const ParticleIndexPair &pi) const {
//   Restraints ret;
//   if (evaluate_index(m, pi, nullptr) < 0) {
//     return Restraints(1, IMP::internal::create_tuple_restraint(this, m, pi));
//   } else {
//     return Restraints();
//   }
// }



IMPNPCTRANSPORT_END_NAMESPACE
