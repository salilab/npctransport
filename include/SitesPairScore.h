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
#include "SitesPairScoreParameters.h"
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
  bool is_orientational_score_; // if true, use orientation-dependent score
  SitesPairScoreParameters params_;
  algebra::Sphere3Ds
    sites0_, sites1_;


  //! Maximal square of distance between particles with interacting sites
  //! based on the interaction range and list of sites
  double ubound_distance2_;

  //! Cache:
  typedef IMP_KERNEL_SMALL_UNORDERED_MAP<ParticleIndex,internal::RigidBodyInfo>
    t_particles_rb_cache;
  mutable t_particles_rb_cache particles_rb_cache_;
  mutable bool is_cache_active_;
  mutable unsigned int cur_cache_id_; // to keep track of caching rounds

 public:

  /**
     For positive sigmas, this is an orientation dependent score between two
     spherical particles that
     contain a fixed set of interaction sites, sites0 and sites1 (for
     first and second particle, resp.).

     The interaction is composed of a non-specific interaction term
     between the bead shperes, and the sum of interactions between
     specific interactions sites of each bead.

     If \sigma{0} and \sigma{1} are positive, the attractive force
     between pairs of sites at an optimal orientation (sites facing
     each other) depends on their distance x. The attraction force
     magnitude is k*x when x<=0.5*range and k*(range-x) when x is
     between 0.5*range and range.  When site0 is rotated by \sigma <
     \sigma{0}, this force decays further by a factor
     (cos{\sigma}-cos{\sigma}0)/(1.0-cos{\sigma0}). The force decays
     similarly when site1 is rotated by \sigma < \sigma{1}.  The
     maximal potential energy difference due to such pair of
     interacting sites is:
       max\DELTA{U}_{site_site} = 0.25 * k * range^2 [kcal/mol]

     If sigma0_deg or sigma1_deg are non-positive, the attractive
     force between pairs of sites within the attraction range is a
     constant k, and the maximal interaction energy is instead:
       max\DELTA{U}_{site-site} = k * range [kCal/mol]

     In addition to site-site interaction, theere is a constanct
     attractive force k_nonspec_attraction between the sphere surfaces
     up to a range range_nonspec_attraction, with maximal energy
     contribution:
       max\DELTA{U}_{non-specifiec} = 0.5 * k_nonspec_attraction * range_nonspec_attraction [kcal/mol]

     Note that for a specific pair of particles, each particle might have
     a different reference frame (rigid body translation and rotation),
     which is applied to the sites list upon score evaluation.

     @param range        Maximal range of site specific attraction in any direction
                         of specific sites placed on particles
     @param k            Maximal site specific attraction coefficient (in
                         kcal/mol/A^2 when sigma0_deg and sigma1_deg are positive,
                         or kcal/mol/A otherwise, i.e., for orientation-independent interacations)
     @param sigma0_deg, sigma1_deg Maximal rotational range of sites 0 and 1, resepctively,
                         on the particle surface, specified in degrees. If either is 0,
			 the pair score between site centers is used with a constant k.
     @param range_nonspec_attraction  Range for non-specific attraction term
                                      between particles that contain the sites
     @param k_nonspec_attraction  Non-specific attraction coefficient between particles
                                  (constant force in kCal/mol/A within specified range)
     @param k_nonspec_repulsion   Repulsion coefficient between particles (constant force
                                  applied when particle spheres overlap in kCal/mol/A)
     @param sites0      List of sites on the first particle, in its local reference frame
     @param sites1      List of sites on the second particle, in its local reference frame
   */
  SitesPairScore(double range, double k,
		 double sigma0_deg, double sigma1_deg,
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

  //  Restraints do_create_current_decomposition(Model *m,
  //                                           const ParticleIndexPair &vt)
  //  const IMP_OVERRIDE;

  //! return the range for site-site attraction
  double get_sites_range() const { return params_.r; }

  //! return the k for site-site attraction
  double get_sites_k() const { return params_.k; }

  SitesPairScoreParameters get_parameters() const {return params_;}

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
