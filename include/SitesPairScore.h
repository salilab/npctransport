/**
 *  \file SitesPairScore.h
 *  \brief A Score on the distance between a pair of particles.
 *
 *  Copyright 2007-2013 IMP Inventors. All rights reserved.
 */

// TODO: verify if energy units are kcal/mol or KT

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
#include <IMP/set_map_macros.h>
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
class IMPNPCTRANSPORTEXPORT SitesPairScore
: public LinearInteractionPairScore
{
 private:
  typedef LinearInteractionPairScore P;


  /************************* class variables ****************/
  bool is_skewed_; // if true, use params_normal_ and params_tangent_
  internal::SitesPairScoreParams
    params_unskewed_; // for any direction (if is_skewed,still gives upper bound)
  internal::SitesPairScoreParams
    params_normal_;  // for normal directiom, if is_skewed_
  internal::SitesPairScoreParams
    params_tangent_; // for tangent directiom, if is_skewed_

  algebra::Sphere3Ds
    sites_,  // sites to be searched against nn_(nnsites_)
    nnsites_;               // sites to be stored in nn_ (nearest neighbours)

  /** sites_first_ is true if sites_ is associated with the first
      particle and nnsites_ with the second particle, false otherwise
      (this switching may be made to improve running time efficiency) */
  bool sites_first_;

  //! Maximal square distance between particles whose sites interact
  //! based on range and list of sites
  double ubound_distance2_;

  //! Cache:
  typedef IMP_KERNEL_SMALL_UNORDERED_MAP<ParticleIndex,internal::RigidBodyInfo>
    t_particles_rb_cache;
  mutable t_particles_rb_cache particles_rb_cache_;
  mutable bool is_cache_active_;
  mutable unsigned int cur_cache_id_; // to keep track of caching rounds

  //! Data structure for finding nearest neighbor (use obsolete)
  //  PointerMember<algebra::NearestNeighbor3D> nn_;

 public:
  /**
     TODO: IMP_DEPRECATED!!!

     A score between two spherical particles that contain a fixed set
     of interaction sites,  sites0 and sites1 (for first and second particle,
     resp.).

     The interaction is decomposed into a non-specific interaction
     between the beads, and the sum of specific interactions between
     sites on each bead.

     The total energetic potential of each individual site-site interaction (in kcal/mol)
       DELTA-U = 0.25 * k * r^2

     This is in addition to contribution from the non-specific interaction (in kcal/mol)
       DELTA-U= 0.5 * k_nonspec_attraction * range_nonspec_attraction

     Note that for a specific pair of particles, each particle might have
     a different reference frame (rigid body translation and rotation),
     which is applied to the sites list upon score evaluation.

     @param range    range of site specific attraction
     @param k        site specific attraction coefficient
     @param range_nonspec_attraction range for non-specific attraction
                                     between particles that contain the sites
     @param k_nonspec_attraction     non-specific attraction coefficient
     @param k_nonspec_repulsion    repulsion coefficient between penetrating particles
     @param sites0    list of sites on the first particle
     @param sites1    list of sites on the second particle

     \deprecated_at{2.2} Use skewed constructor instead, with 1.0 values if
     unskewedness needed.
*/
  IMPKERNEL_DEPRECATED_METHOD_DECL(2.2)
  SitesPairScore(double range, double k,
                 double range_nonspec_attraction,
                 double k_nonspec_attraction,
                 double k_nonspec_repulsion,
                 const algebra::Sphere3Ds &sites0,
                 const algebra::Sphere3Ds &sites1);

  /**
     A skered version for a score between two spherical particles that
     contain a fixed set of interaction sites, sites0 and sites1 (for
     first and second particle, resp.).

     The interaction is decomposed into a non-specific interaction
     between the beads, and the sum of specific interactions between
     sites on each bead.

     The skewing is with regard to directionality of the site-site
     interaction. The interaction vector between a pair of sites is decomposed into a
     surface-normal component and a tangent component, each with a
     different range and force constant, as specified by range_skew
     and k_skew. The surface normal is the vector between the two bead centers.

     The total energetic potential for each pair of site-site interactions, regardless of skewing, is (in kcal/mol):
       DELTA-U = 0.25 * k * r^2

     This is in addition to contribution from the non-specific interaction (in kcal/mol)
       DELTA-U= 0.5 * k_nonspec_attraction * range_nonspec_attraction

     Note that for a specific pair of particles, each particle might have
     a different reference frame (rigid body translation and rotation),
     which is applied to the sites list upon score evaluation.

     @param range       Maximal range of site specific attraction in any direction
                         of specific sites placed on particles
     @param k           Site specific attraction coefficient
                         (combination of normal and tangent contributions)
     @param range2_skew  r_tangent^2/r_normal^2 ratio  (s.t. r_tangent^2 + r_normal^2 = sites_range^2)
     @param k_skew      k_tangent/k_normal ratio  (s.t. k_tangent * k_normal = 4*k_attraction)
     @param range_nonspec_attraction  Range for non-specific attraction
                                      between particles that contain the sites
     @param k_nonspec_attraction  Non-specific attraction coefficient between particles
     @param k_nonspec_repulsion   Repulsion coefficient between penetrating particles
     @param sites0      List of sites on the first particle
     @param sites1      List of sites on the second particle
   */
  SitesPairScore(double range, double k,
                 double range2_skew, double k_skew,
                 double range_nonspec_attraction, double k_nonspec_attraction,
                 double k_nonspec_repulsion,
                 const algebra::Sphere3Ds &sites0,
                 const algebra::Sphere3Ds &sites1);


 public:

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
    ( Model *m, const ParticleIndexPairs &p,
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
  double get_sites_range() const { return params_unskewed_.r; }

  //! return the k for site-site attraction
  double get_sites_k() const { return params_unskewed_.k; }


 public:
  IMP_OBJECT_METHODS(SitesPairScore);

 private:
  // gets the rigid body information (e.g., translation, inverse rotation)
  // associated with particle m.pi, possibly from cache (depending on internal
  // cache definitions)
  inline internal::RigidBodyInfo
    get_rigid_body_info(Model* m, ParticleIndex pi) const;

 private:
  // sets the sites associated with each partner to sites0
  // and sites1, respectively (in local reference frame)
  void set_sites(const algebra::Sphere3Ds &sites0,
                 const algebra::Sphere3Ds &sites1);

 private:
  // maintain a cache for evaluate_index() till call to deactivate_cache()
  inline void activate_cache() const
  {
    // (note: cache is mutable)
    is_cache_active_ = true;  cur_cache_id_++;
  }

  inline void deactivate_cache() const
  {
    // (note: cache is mutable)
    is_cache_active_ = false;
  }


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



IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_SITES_PAIR_SCORE_H */
