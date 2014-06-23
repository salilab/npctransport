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

IMPNPCTRANSPORT_BEGIN_NAMESPACE

namespace {
  // return maximum sum of radii of pairs of vectors in R0 and R1
  double get_max_r_sum(const algebra::Vector3Ds& R0,
                const algebra::Vector3Ds& R1) {
    double max = 0.0;
    for(unsigned int i = 0; i < R0.size(); i++){
      for(unsigned int j = 0; j < R1.size(); j++){
        double s = R0[i].get_magnitude() + R1[j].get_magnitude();
        max = (max > s) ? max : s;
      }
    }
    return max;
  }
}

/**
   A score between two spherical particles that contain a fixed set
   of interaction sites,  sites0 and sites1 (for first and second particle,
   resp.).

   Note that for a specific pair of particles, each particle might have
   a different reference frame (rigid body translation and rotation),
   which is applied to the sites list upon score evaluation.

   @param range          range of site specific attraction
   @param k_attraction   site specific attraction coefficient
   @param range_nonspec_attraction range for non-specific attraction
                                   between particles that contain the sites
   @param k_nonspec_attraction     non-specific attraction coefficient
   @param k_repulsion    repulsion coefficient between penetrating particles
   @param sites0         list of sites on the first particle
   @param sites1         list of sites on the second particle
*/
SitesPairScore::SitesPairScore(double sites_range, double k_attraction,
                               double range_nonspec_attraction,
                               double k_nonspec_attraction, double k_repulsion,
                               const algebra::Vector3Ds &sites0,
                               const algebra::Vector3Ds &sites1)
    : P(k_repulsion, range_nonspec_attraction, k_nonspec_attraction,
        "SitesPairScore %1%"),
      sites_range_(sites_range),
      sites_k_(k_attraction),
      is_cache_active_(false),
      cur_cache_id_(INVALID_CACHE_ID)
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
  nn_ = new algebra::NearestNeighbor3D(nnsites_);
  max_sites_r_sum_ = get_max_r_sum(sites0, sites1);
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
  LinearInteractionPairScore::EvaluationCache lips_cache =
    P::get_evaluation_cache();

  // II. Return if parent particles are out of site-specific interaction range
  //     using cache to avoid some redundant calcs
  double particles_delta = std::sqrt(lips_cache.particles_delta_squared);
  double min_sites_delta = particles_delta - max_sites_r_sum_;
  IMP_LOG(TERSE, "min_sites_delta " << min_sites_delta
          << " ; range " << sites_range_ << std::endl);  // TODO: VERBOSE
  if (min_sites_delta > sites_range_) {
    IMP_LOG(TERSE, "Sites contribution is 0.0 and non-specific score  is "
                       << non_specific_score << std::endl);  // TODO: VERBOSE
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
  double const sites_range2 = algebra::get_squared(sites_range_);
  for (unsigned int i = 0; i < sites_.size(); ++i) {
    // filter to evaluate only sites within range of attraction:
    // algebra::Vector3D trp = relative.get_transformed(sites_[i]);
    //    Ints nn = nn_->get_in_ball(trp, sites_range_);
    //for (unsigned int j = 0; j < nn.size(); ++j) {
    algebra::Vector3D g0 = rbi0.tr.get_transformed(sites_[i]);
    IMP_LOG_PROGRESS( "g0 = " << g0 );
    for(unsigned int j = 0 ; j < nnsites_.size(); ++j) {
      // double d2=algebra::get_squared(range_);
      algebra::Vector3D g1 = rbi1.tr.get_transformed(nnsites_[j]);
      IMP_LOG_PROGRESS( "Evaluating sites " << g0 << " ; " << g1 << std::endl );
      sum +=
        internal::evaluate_one_site_3(sites_k_,
                                      sites_range_, sites_range2,
                                      rbi0, rbi1,
                                      sites_[i], nnsites_[j],
                                      g0, g1,
                                      da);
      IMP_LOG_PROGRESS( "Sum " << sum << std::endl);
    }
  }
  IMP_LOG(VERBOSE,
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


#define IMP_HS(na, nb)                                   \
  else if (sites0.size() == na && sites1.size() == nb) { \
    typedef TemplateSitesPairScore<na, nb, true> TSPS;   \
    IMP_NEW(TSPS, ps, (range, ka, kr, sites0, sites1));  \
    ppr->set_score(value, ps.get());                     \
  }

#if IMP_BUILD == IMP_FAST
#define IMP_S0 \
  (1)(2)(3)(4)(5)(6)(7)(8)(9)(10)(11)(12)(13)(14)(15)(16)(17)(18)(19)(20)
#else
#define IMP_S0 (1)(2)(3)(4)(5)
#endif
#define IMP_CALL(a, bb) a(BOOST_PP_SEQ_ELEM(0, bb), BOOST_PP_SEQ_ELEM(1, bb))

#define IMP_MACRO(r, product) IMP_CALL(IMP_HS, product)

/*
// #define IMP_ADD_CASE(na,nb)  else if (sites0.size()==na
//                                       && sites1.size()==nb) {
//     typedef TemplateSitesPairScore<na, nb, true> TSPS;
//     IMP_NEW(TSPS, ps, (range, ka, rangena, kna, kr, sites0, sites1));
//     ppr->set_score(value, ps.get());
//   }
*/

/**
   @param rangesa  range of specific attraction
   @param ksa      coefficient for specific attraction
   @param rangensa range for non-specific attraction
   @param knsa      coefficient for non-specific attraction
   @param kr      coefficient for repulsion between penetrating particles
   @param sites0   sites on side of first particle
   @param sites1   sites on side of other particle
*/
IMP::PairScore* create_sites_pair_score
( double rangesa, double ksa, double rangensa,
  double knsa, double kr,
  const algebra::Vector3Ds &sites0,
  const algebra::Vector3Ds &sites1)
{
  // use the point-location based one if appropriate
  if ((sites0.size() >= 4 && sites1.size() >= 4 &&
       (sites0.size() > .5 * sites1.size() ||
        sites1.size() > .5 * sites0.size()))) {
    IMP_NEW(SitesPairScore, ps,
            (rangesa, ksa, rangensa, knsa, kr, sites0, sites1));
    return ps.release();
  }
  // TODO: check if the approach of the next clause (originally
  // IMP_ADD_CASE) can accelerate computations by having a template
  // for spefified site sizes
  // BOOST_PP_SEQ_FOR_EACH_PRODUCT(IMP_MACRO, (IMP_S0)(IMP_S0))
  //  IMP_ADD_CASE(1,15)
  else if (sites0.size() == 1 && sites1.size() == 15) {
    typedef TemplateSitesPairScore<1, 15, true> TSPS;
    IMP_NEW(TSPS, tsps, (rangesa, ksa, rangensa, knsa, kr,
                         sites0, sites1));
    return tsps.release();
  } else {
    IMP_NEW(SitesPairScore, ps, (rangesa, ksa, rangensa, knsa, kr,
                                 sites0, sites1));
    return ps.release();
  }
}

IMPNPCTRANSPORT_END_NAMESPACE
