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
        max = (max > s) ? max : s;
      }
    }
    return max;
  }
}

// unskewed version of interaction (indifferent to interaction directionality)
SitesPairScore::SitesPairScore(double range, double k,
                               double range_nonspec_attraction,
                               double k_nonspec_attraction,
                               double k_nonspec_repulsion,
                               const algebra::Sphere3Ds &sites0,
                               const algebra::Sphere3Ds &sites1)
    : P(k_nonspec_repulsion, range_nonspec_attraction,
        k_nonspec_attraction, "SitesPairScore %1%"),
      is_skewed_(false),
      params_unskewed_(range, k),
      params_normal_(0.0, 0.0),
      params_tangent_(0.0, 0.0),
      is_cache_active_(false),
      cur_cache_id_(INVALID_CACHE_ID)
{
  set_sites(sites0, sites1);
}

// skewed version of interaction score (different in
// normal and tangent directions)
SitesPairScore::SitesPairScore(double range, double k,
                               double range2_skew, double k_skew,
                               double range_nonspec_attraction,
                               double k_nonspec_attraction,
                               double k_repulsion,
                               const algebra::Sphere3Ds &sites0,
                               const algebra::Sphere3Ds &sites1)
  :
  P(k_repulsion, range_nonspec_attraction, k_nonspec_attraction,
    "SitesPairScore %1%"),
  params_unskewed_(range, k),
  params_normal_(0.0, 0.0),
  params_tangent_(0.0, 0.0),
  is_cache_active_(false),
  cur_cache_id_(INVALID_CACHE_ID)
{
  // IMP_ALWAYS_CHECK(range2_skew > 0.0 && k_skew > 0.0,
  //                  "Skew must be positive", base::ValueException);
  is_skewed_ = (range2_skew != 0.0 && k_skew != 0.0);
  if(is_skewed_){
    double kN = sqrt(k / k_skew);
    double kT = sqrt(k * k_skew);
    double const r2 = range*range;
    double const& s = range2_skew; // shorthand
    double rN = sqrt(r2/(s+1));
    double rT = sqrt(s*r2/(s+1));
    params_normal_.set_rk(rN, kN);
    params_tangent_.set_rk(rT, kT);
    IMP_LOG_PROGRESS("Ks: " << k
                     << " ; kN " << kN << " , " << params_normal_.k
                     << " ; kT " << kT << " , " << params_tangent_.k
                     << std::endl);
    IMP_LOG_PROGRESS("Range: " << range
                     << " ; rN " << rN << " , " << params_normal_.r
                     << " ; rT " << rT << " , " << params_tangent_.r
                     << std::endl);
  } else {
    std::cout << "Creating old version of SitesPairScore, for backward compatibility" << std::endl;
  }
  set_sites(sites0, sites1);
}


void SitesPairScore::set_sites(const algebra::Sphere3Ds &sites0,
               const algebra::Sphere3Ds &sites1)
{
  // store the big set of sizes in nnsites_
  if (sites0.size() > sites1.size()) {
    sites_first_ = false;
    sites_ = sites1;
    nnsites_ = sites0;
  } else {
    sites_first_ = true;
    sites_ = sites0;
    nnsites_ = sites1;
  }
  //  nn_ = new algebra::NearestNeighbor3D(nnsites_); // Need to convert to get_center() of each
  // Find upper bound for distance between particles whose sites interact
  double ubound_distance = (get_max_r_sum(sites0, sites1) + params_unskewed_.r);
  this->ubound_distance2_ = ubound_distance * ubound_distance;
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
  // make sure rb0 is associated with sites_, and rb1 with nn_sites_:
  if (!sites_first_) std::swap(pi0, pi1);
  // get rbi0/1 from cache, update if needed
  internal::RigidBodyInfo rbi0 = get_rigid_body_info(m, pi0);
  internal::RigidBodyInfo rbi1 = get_rigid_body_info(m, pi1);
  IMP_LOG_PROGRESS( "RBI0.rb " << rbi0.rb
                    << " RB0.cache_id " << rbi0.cache_id
                    << "RBI0.tr " << rbi0.tr << std::endl);
  IMP_LOG_PROGRESS( "RBI1.cache_id " << rbi1.cache_id
                    << "RBI1.RB " << rbi1.rb
                    << "RBI1.tr " << rbi1.tr << std::endl);
  // bring sites_ to the correct orientation relative to nn_sites_[j]
  // algebra::Transformation3D relative = c1->tr.get_inverse() * c0->tr;
  // sum over specific interactions between all pairs of sites:
  double sum = 0;
  for (unsigned int i = 0; i < sites_.size(); ++i) {
    // filter to evaluate only sites within range of attraction:
    // algebra::Vector3D trp = relative.get_transformed(sites_[i]);
    //    Ints nn = nn_->get_in_ball(trp, sites_range_);
    //for (unsigned int j = 0; j < nn.size(); ++j) {
    algebra::Vector3D g0 = rbi0.tr.get_transformed(sites_[i].get_center());
    IMP_LOG_PROGRESS( "g0 = " << g0 );
    for(unsigned int j = 0 ; j < nnsites_.size(); ++j) {
      // double d2=algebra::get_squared(range_);
      algebra::Vector3D g1 = rbi1.tr.get_transformed(nnsites_[j].get_center());
      IMP_LOG_PROGRESS( "Evaluating sites " << g0 << " ; " << g1 << std::endl );
      if(is_skewed_)
        {
          sum +=
            internal::evaluate_one_site_4(params_unskewed_,
                                          params_normal_,
                                          params_tangent_,
                                          rbi0, rbi1,
                                          sites_[i], nnsites_[j],
                                          g0, g1,
                                          da);
      } else
        {
          sum +=
            internal::evaluate_one_site_3(params_unskewed_.k,
                                          params_unskewed_.r,
                                          rbi0, rbi1,
                                          sites_[i], nnsites_[j],
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

Restraints SitesPairScore::do_create_current_decomposition(
    Model *m, const ParticleIndexPair &pi) const {
  Restraints ret;
  if (evaluate_index(m, pi, nullptr) < 0) {
    return Restraints(1, IMP::internal::create_tuple_restraint(this, m, pi));
  } else {
    return Restraints();
  }
}



IMPNPCTRANSPORT_END_NAMESPACE
