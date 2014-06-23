/**
 *  \file SitesPairScore.h
 *  \brief A Score on the distance between a pair of particles.
 *
 *  Copyright 2007-2013 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_SITES_PAIR_SCORE_H
#define IMPNPCTRANSPORT_SITES_PAIR_SCORE_H

#include "npctransport_config.h"
#include "linear_distance_pair_scores.h"
#include "internal/RigidBodyInfo.h"
#include "internal/sites.h"
#include <IMP/PairScore.h>
#include <IMP/UnaryFunction.h>
#include <IMP/Pointer.h>
#include <IMP/core/XYZR.h>
#include <IMP/core/Typed.h>
#include <IMP/core/SphereDistancePairScore.h>
#include <IMP/display/particle_geometry.h>
#include <IMP/generic.h>
#include <IMP/algebra/vector_search.h>
#include <IMP/algebra/Transformation3D.h>
#include <IMP/base/set_map_macros.h>
#include <IMP/container/PredicatePairsRestraint.h>
#include <IMP/atom/estimates.h>
#include <boost/unordered_set.hpp>

#include <boost/array.hpp>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

/** \brief  Apply a function to the distance between two particles with
            a set of specific binding sites

    The sites are expressed in the local reference frame of
    the two rigid bodies. Care must be taken to pass the bodies
    in the appropriate order. See construction documentation for more details.
*/
class IMPNPCTRANSPORTEXPORT SitesPairScore : public LinearInteractionPairScore {
  typedef LinearInteractionPairScore P;

  // range and coefficient of site specific attraction:
  double sites_range_, sites_k_;

  algebra::Vector3Ds
    sites_,  // sites to be searched against nn_(nnsites_)
    nnsites_;               // sites to be stored in nn_ (nearest neighbours)

  // true if sites_ is associated with the first particle and nnsites_ with the
  // second particle, false otherwise
  // (this switching may be made to improve running time efficiency)
  bool sites_first_;

  // the maximal sum of radii of any pair of sites
  // from sites_ and nnsites_ respectively
  double max_sites_r_sum_;

 private:
  typedef IMP_BASE_SMALL_UNORDERED_MAP<ParticleIndex,internal::RigidBodyInfo>
    t_particles_rb_cache;
  mutable t_particles_rb_cache particles_rb_cache_;
  mutable bool is_cache_active_;
  mutable unsigned int cur_cache_id_; // to keep track of caching rounds

  PointerMember<algebra::NearestNeighbor3D> nn_;

 private:
  // maintain a cache for evaluate_index() till
  // call to deactivate_cache()
  void activate_cache() const // cache is mutable
  {
    is_cache_active_ = true;
    cur_cache_id_++; // invalidate old entries
    //    std::cout << "Activating cache with cache id " << cur_cache_id_ << std::endl;
  }

  void deactivate_cache() const // cache is mutable
  {
    is_cache_active_ = false;
    //    std::cout << "Deactivating cache" << std::endl;
  }

 public:
  /**
     A score between two spherical particles that contain a fixed set
     of interaction sites,  sites0 and sites1 (for first and second particle,
     resp.).

     Note that for a specific pair of particles, each particle might have
     a different reference frame (rigid body translation and rotation),
     which is applied to the sites list upon score evaluation.

     @param sites_range    range of site specific attraction
     @param k_attraction   site specific attraction coefficient
     @param range_nonspec_attraction range for non-specific attraction
                                     between particles that contain the sites
     @param k_nonspec_attraction     non-specific attraction coefficient
     @param k_repulsion    repulsion coefficient between penetrating particles
     @param sites0         list of sites on the first particle
     @param sites1         list of sites on the second particle
   */
  SitesPairScore(double sites_range, double k_attraction,
                 double range_nonspec_attraction, double k_nonspec_attraction,
                 double k_repulsion, const algebra::Vector3Ds &sites0,
                 const algebra::Vector3Ds &sites1);

  virtual double evaluate_indexes(Model *m, const ParticleIndexPairs &p,
                                  DerivativeAccumulator *da,
                                  unsigned int lower_bound,
                                  unsigned int upper_bound) const IMP_FINAL {
    activate_cache();
    double ret=0.0;
    for (unsigned int i = lower_bound; i < upper_bound; ++i) {
      ret += evaluate_index(m, p[i], da);
    }
    deactivate_cache();
    return ret;
  }

  //! evaluated indexes for the range from lower_bound to upper_bound
  //! in p, if score>max then return max value of double
  double evaluate_if_good_indexes
    ( kernel::Model *m, const kernel::ParticleIndexPairs &p,
      DerivativeAccumulator *da,
      double max, unsigned int lower_bound, unsigned int upper_bound) const
  {
    activate_cache();
    double ret = 0.0;
    for (unsigned int i = lower_bound; i < upper_bound; ++i) {
      ret += evaluate_if_good_index(m, p[i], da, max - ret);
      if (ret > max) return std::numeric_limits<double>::max();
    }
    deactivate_cache();
    return ret;
  }


  /** evaluate the score for the pair of model particle indexes in p,
      updating score derivatives to da
  */
  virtual double evaluate_index(Model *m, const ParticleIndexPair &p,
                                DerivativeAccumulator *da) const IMP_OVERRIDE;

  virtual ModelObjectsTemp do_get_inputs(Model *m,
                                         const ParticleIndexes &pis) const;

  Restraints do_create_current_decomposition(Model *m,
                                             const ParticleIndexPair &vt)
    const IMP_OVERRIDE;

  //! return the range for site-site attraction
  double get_sites_range() const { return sites_range_; }

  //! return the k for site-site attraction
  double get_sites_k() const { return sites_k_; }


 public:

  IMP_OBJECT_METHODS(SitesPairScore);

 private:

  // gets the rigid body information (e.g., translation, inverse rotation)
  // associated with particle m.pi, possibly from cache (depending on internal
  // cache definitions)
  inline internal::RigidBodyInfo
    get_rigid_body_info(Model* m, ParticleIndex pi) const;


};

internal::RigidBodyInfo
SitesPairScore::get_rigid_body_info
(Model* m, ParticleIndex pi) const
{
  IMP_USAGE_CHECK(core::RigidBody::get_is_setup(m, pi),
                  "PI " << pi.get_index() << " not a rigid body");
  if(is_cache_active_){
    std::pair<t_particles_rb_cache::iterator, bool>
      p = particles_rb_cache_.insert
      (std::make_pair(pi, internal::RigidBodyInfo()));
    internal::RigidBodyInfo& rbi_cached = p.first->second;
    bool const rbi_was_in_cache = !p.second;
    if(!rbi_was_in_cache || rbi_cached.cache_id != cur_cache_id_)
      {
        rbi_cached.set_particle(m, pi, cur_cache_id_);
      }
    return rbi_cached;
  } else {
    return internal::RigidBodyInfo(m, pi, INVALID_CACHE_ID);
  }
}


class TemplateBaseSitesPairScore : public LinearInteractionPairScore {
 protected:
  typedef LinearInteractionPairScore P;
  double range_, k_;
  inline static boost::tuple<algebra::Transformation3D,
                             algebra::Transformation3D, algebra::Rotation3D,
                             algebra::Rotation3D>
  get_transformations(core::RigidBody rb0, core::RigidBody rb1) {
    algebra::Transformation3D tr0 =
        rb0.get_reference_frame().get_transformation_to();
    algebra::Transformation3D tr1 =
        rb1.get_reference_frame().get_transformation_to();
    algebra::Rotation3D itr0 = tr0.get_rotation().get_inverse();
    algebra::Rotation3D itr1 = tr1.get_rotation().get_inverse();
    // algebra::Transformation3D relative=tr1.get_inverse() *tr0;
    return boost::make_tuple(tr0, tr1, itr0, itr1);
  }

 public:
  TemplateBaseSitesPairScore(double range, double k_attraction,
                             double range_nonspec_attraction,
                             double k_nonspec_attraction, double k_repulsion,
                             std::string name)
      : P(k_repulsion, range_nonspec_attraction, k_nonspec_attraction, name),
        range_(range),
        k_(k_attraction) {}
  //! return the upper bound on the range over which it is non-zero
  Restraints do_create_current_decomposition(Model *m,
                                             const ParticleIndexPair &vt) const
      IMP_OVERRIDE {
    Restraints ret;
    if (evaluate_index(m, vt, nullptr) < 0) {
      return Restraints(1, IMP::internal::create_tuple_restraint(this, m, vt));
    } else {
      return Restraints();
    }
  }
};

/**
   pair score with NA sites vs. NB sites in each pair

   WHICH is a variant for which type of evaluation is used (advanced,
   use True if you don't know better)
 */
template <unsigned int NA, unsigned int NB, bool WHICH>
class TemplateSitesPairScore : public TemplateBaseSitesPairScore {
  typedef TemplateBaseSitesPairScore P;
  boost::array<algebra::Vector3D, NA> sites0_;
  boost::array<algebra::Vector3D, NB> sites1_;

 public:
  TemplateSitesPairScore(double range, double k_attraction,
                         double range_nonspec_attraction,
                         double k_nonspec_attraction, double k_repulsion,
                         const algebra::Vector3Ds &sites0,
                         const algebra::Vector3Ds &sites1)
      : P(range, k_attraction, k_repulsion, range_nonspec_attraction,
          k_nonspec_attraction, "Sites %1%") {
    std::copy(sites0.begin(), sites0.end(), &sites0_[0]);
    std::copy(sites1.begin(), sites1.end(), &sites1_[0]);
  }
  inline double evaluate_index(Model *m, const ParticleIndexPair &pp,
                               DerivativeAccumulator *da) const {
    core::RigidBody rb0(m, pp[0]), rb1(m, pp[1]);
    algebra::Transformation3D tr0, tr1;
    algebra::Rotation3D itr0, itr1;
    boost::tie(tr0, tr1, itr0, itr1) = P::get_transformations(rb0, rb1);
    double sum = 0;
    for (unsigned int i = 0; i < NA; ++i) {
      for (unsigned int j = 0; j < NB; ++j) {
        // double d2=algebra::get_squared(range_);
        if (WHICH) {
          sum += internal::evaluate_one_site(P::k_, P::range_, rb0, rb1, tr0,
                                             tr1, itr0, itr1, sites0_[i],
                                             sites1_[j], da);
        } else {
          sum += internal::evaluate_one_site_2(P::k_, P::range_, rb0, rb1, tr0,
                                               tr1, itr0, itr1, sites0_[i],
                                               sites1_[j], da);
        }
      }
    }
    return sum + P::evaluate_index(m, pp, da);
  }
};

/**
   Create a score function between particles with site lists
   sites0 and sites1 on their surface, with specified interaction
   parameters.

   @param rangesa  range for specific attraction
   @param ksa      coefficient for specific attraction
   @param rangensa range for non-specific attraction
   @param knsa     coefficient for non-specific attraction
   @param kr       coefficient for repulsion between penetrating particles
   @param sites0   location of sites on particles from one side
   @param sites1   location of sites on particles from other side

   @return a PairScore that works on site-site interaction between
           adjacent particles
*/
IMPNPCTRANSPORTEXPORT
IMP::PairScore* create_sites_pair_score
( double rangesa, double ksa, double rangensa, double knsa,
  double kr, const algebra::Vector3Ds &sites0,
  const algebra::Vector3Ds &sites1);


IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_SITES_PAIR_SCORE_H */
