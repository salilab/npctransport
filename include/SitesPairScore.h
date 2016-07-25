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
  typedef IMP_KERNEL_LARGE_UNORDERED_MAP<ParticleIndex,internal::RigidBodyInfo>
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
                                  unsigned int upper_bound) const IMP_FINAL;

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

#ifndef SWIG
  /**
     EvaluatE all site-site interactions
     for evaluate_index() for the pair pip in model m. If da is not nullptr,
     it accumulated appropriate derivatives. If n_contacts_accumulator>=0,
     then the number of individual contacts found is asscumulated there.

     @param m the model
     @param pip the pair of particle indexes in m
     @param da optional accumulator for force and torque derivatives
     @param n_contacts_accumulator optional pointer to accumulator for number of
           individual contacts.

     @return the site-site contributions for the score for the pair
             pip in model m.
  */
  double
    evaluate_site_contributions_with_internal_tables
    (algebra::Sphere3D const* spheres_table,
     double const**quaternions_tables,
     algebra::Sphere3D *sphere_derivatives_table,
     double **torques_tables,
     const ParticleIndexPair &pip,
     DerivativeAccumulator *da,
     int* n_contacts_accumulator = nullptr
     ) const;

  /**
     EvaluatE all site-site interactions
     for evaluate_index() for the pair pip in model m. If da is not nullptr,
     it accumulated appropriate derivatives. If n_contacts_accumulator>=0,
     then the number of individual contacts found is asscumulated there.

     @param m the model
     @param pip the pair of particle indexes in m
     @param da optional accumulator for force and torque derivatives
     @param n_contacts_accumulator optional pointer to accumulator for number of
           individual contacts.

     @return the site-site contributions for the score for the pair
             pip in model m.
  */
  double evaluate_site_contributions
    (Model* m,
     const ParticleIndexPair &pip,
     DerivativeAccumulator *da,
     int* n_contacts_accumulator = nullptr) const;

#endif

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

  /** evaluate the score for the pair of model particle indexes in p,
      updating score derivatives to da, and using internal attribute
      tables in Model
  */
  inline double evaluate_index_with_internal_tables
    ( Model* m,
      algebra::Sphere3D const* spheres_table,
      double const **quaternions_tables,
      algebra::Sphere3D *sphere_derivatives_table,
      double **torques_tables,
      const ParticleIndexPair &p,
      DerivativeAccumulator *da) const;

  // gets the rigid body information (e.g., translation, inverse rotation)
  // associated with particle m.pi, possibly from cache (depending on internal
  // cache definitions)
  inline internal::RigidBodyInfo
    get_rigid_body_info
    (algebra::Sphere3D const* spheres_table,
     double const** quaternions_tables,
     ParticleIndex pi) const;

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

//!
inline double
SitesPairScore::evaluate_index
(Model *m, const ParticleIndexPair &p,
 DerivativeAccumulator *da) const{
  // get internal tables:
  algebra::Sphere3D const* spheres_table=
    m->access_spheres_data();
  double const* quaternions_tables[4];
  for(unsigned int i = 0; i < 4; i++){
    quaternions_tables[i]=
      core::RigidBody::access_quaternion_i_data(m, i);
  }
  algebra::Sphere3D* sphere_derivatives_table=
    m->access_sphere_derivatives_data();
  double* torques_tables[3];
  for(unsigned int i = 0; i < 3; i++){
    torques_tables[i]=
      core::RigidBody::access_torque_i_data(m, i);
  }
  // evaluate:
  return evaluate_index_with_internal_tables(m,
                                             spheres_table,
                                             quaternions_tables,
                                             sphere_derivatives_table,
                                             torques_tables,
                                             p,
                                             da);
}




/**
   the sites of each particle are transformed to a common frame of reference
   (using the reference frame of each particle), and the site-specific
   attraction, and the inter-particle non specific attraction and repulsion
   are evaluated and summed.
*/
inline double
SitesPairScore::evaluate_index_with_internal_tables
( Model* m,
  algebra::Sphere3D const* spheres_table,
 double const** quaternions_tables,
 algebra::Sphere3D *sphere_derivatives_table,
 double **torques_tables,
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

  double site_score=evaluate_site_contributions_with_internal_tables
    (spheres_table,
     quaternions_tables,
     sphere_derivatives_table,
     torques_tables,
     pip, da);
  // III. evaluate site-specific contributions :
  return site_score + non_specific_score;
}

#ifndef SWIG

//!
inline double
SitesPairScore::evaluate_site_contributions_with_internal_tables
( algebra::Sphere3D const* spheres_table,
  double const**quaternions_tables,
  algebra::Sphere3D *sphere_derivatives_table,
  double **torques_tables,
  const ParticleIndexPair &pip,
  DerivativeAccumulator *da,
  int* n_contacts_accumulator
  ) const
{
  IMP_OBJECT_LOG;

  // bring sites_ to the frame of reference of nn_sites_ and nn_
  ParticleIndex pi0 = pip[0];
  ParticleIndex pi1 = pip[1];
  // get rbi0/1 from cache, update if needed
  internal::RigidBodyInfo rbi0 = get_rigid_body_info(spheres_table,
                                                     quaternions_tables,
                                                     pi0);
  internal::RigidBodyInfo rbi1 = get_rigid_body_info(spheres_table,
                                                     quaternions_tables,
                                                     pi1);
  IMP_LOG_PROGRESS( "RBI0.pi " << rbi0.pi
                    << " RB0.cache_id " << rbi0.cache_id
                    << "RBI0.tr " << rbi0.tr << std::endl);
  IMP_LOG_PROGRESS( "RBI1.cache_id " << rbi1.cache_id
                    << "RBI1.pi " << rbi1.pi
                    << "RBI1.tr " << rbi1.tr << std::endl);
  // sum over specific interactions between all pairs of sites:
  double sum = 0;
  if(is_orientational_score_){
    // Pre-compute a few variables that do not depend on either both sites or on site1
    algebra::Vector3D const& gRB0= rbi0.tr.get_translation();
    algebra::Vector3D const& gRB1= rbi1.tr.get_translation();
    algebra::Vector3D gUnitRB0RB1= gRB1-gRB0;
    double distRB0RB1= get_magnitude_and_normalize_in_place(gUnitRB0RB1); // distance between centers
    for (unsigned int i = 0; i < sites0_.size(); ++i) {
      algebra::Vector3D gSite0 = rbi0.tr.get_transformed(sites0_[i].get_center());
      algebra::Vector3D gUnitRB0Site0= (gSite0-gRB0)*rbi0.iradius;
      double cosSigma0 = gUnitRB0Site0*gUnitRB0RB1;
      if(cosSigma0 < params_.cosSigma1_max) { // not in range... - note the indexing is not an error - sigma0 is equivalent to params_.sigma1
        continue;
      }
      double kFactor0=internal::get_k_factor(cosSigma0, params_.cosSigma1_max); // note the indexing is not an error - sigma0 is equivalent to params_.sigma1
      algebra::Vector3D gRotSigma0;
      double dKFactor0;
      if(da){
        gRotSigma0 = get_vector_product(gUnitRB0Site0,gUnitRB0RB1);
        double absSinSigma0 = get_magnitude_and_normalize_in_place(gRotSigma0);
        dKFactor0=internal::get_derivative_k_factor(absSinSigma0, params_.cosSigma1_max);
      }
      for(unsigned int j = 0 ; j < sites1_.size(); ++j) {
        algebra::Vector3D gSite1 = rbi1.tr.get_transformed(sites1_[j].get_center());
        IMP_LOG_PROGRESS( "Evaluating sites at global coordinates: " << gSite0
                          << " ; " << gSite1 << std::endl );
        double cur_score;
        cur_score =
          internal::evaluate_pair_of_sites(params_,
                                           rbi0, rbi1,
                                           gSite1,
                                           gUnitRB0RB1, distRB0RB1,
                                           gRotSigma0,
                                           kFactor0, dKFactor0,
                                           da,
                                           sphere_derivatives_table,
                                           torques_tables);
        sum += cur_score;
        if(n_contacts_accumulator){
          (*n_contacts_accumulator) += (cur_score!=0.0);
        }
      } // j
    } // i
  } // is_orientational_score_
  else
    {
      for (unsigned int i = 0; i < sites0_.size(); ++i) {
        algebra::Vector3D g0 = rbi0.tr.get_transformed(sites0_[i].get_center());
        for(unsigned int j = 0 ; j < sites1_.size(); ++j) {
          algebra::Vector3D g1 = rbi1.tr.get_transformed(sites1_[j].get_center());
          IMP_LOG_PROGRESS( "Evaluating sites at global coordinates: " << g0 << " ; " << g1 << std::endl );
          double cur_score;
          // old score
          cur_score =
            internal::evaluate_one_site_3(params_.k,
                                          params_.r,
                                          rbi0, rbi1,
                                          sites0_[i], sites1_[j],
                                          g0, g1,
                                          da,
                                          sphere_derivatives_table,
                                          torques_tables);

          sum += cur_score;
          if(n_contacts_accumulator){
            (*n_contacts_accumulator) += (cur_score!=0.0);
          }
        }// j
      }// i
    }
  IMP_LOG_PROGRESS( "Sum " << sum << std::endl);
  return sum;
}


//!
inline double
SitesPairScore::evaluate_site_contributions
(Model* m,
 const ParticleIndexPair &pip,
 DerivativeAccumulator *da,
 int* n_contacts_accumulator) const
{
  // Get internal tables
  algebra::Sphere3D const* spheres_table=
    m->access_spheres_data();
  double const* quaternions_tables[4];
  for(unsigned int i = 0; i < 4; i++){
    quaternions_tables[i]=
      core::RigidBody::access_quaternion_i_data(m, i);
  }
  algebra::Sphere3D* sphere_derivatives_table=
    m->access_sphere_derivatives_data();
  double* torques_tables[3];
  for(unsigned int i = 0; i < 3; i++){
    torques_tables[i]=
      core::RigidBody::access_torque_i_data(m, i);
  }
  // evaluate:
  return evaluate_site_contributions_with_internal_tables
    (spheres_table,
     quaternions_tables,
     sphere_derivatives_table,
     torques_tables,
     pip,
     da,
     n_contacts_accumulator);
}

#endif // ifndef SWIG

//!
inline internal::RigidBodyInfo
SitesPairScore::get_rigid_body_info
(algebra::Sphere3D const* spheres_table,
 double const **quaternions_tables,
 ParticleIndex pi) const
{
  // TODO: add usage check that it has valid quaternions
  //  IMP_USAGE_CHECK(core::RigidBody::get_is_setup(m, pi),
  //                "PI " << pi.get_index() << " not a rigid body");
  if(is_cache_active_){
    std::pair<t_particles_rb_cache::iterator, bool>
      p = particles_rb_cache_.insert
      (std::make_pair(pi, internal::RigidBodyInfo()));
    internal::RigidBodyInfo& rbi_cached = p.first->second;
    bool const rbi_in_cache = !p.second;
    if(!rbi_in_cache ||
       rbi_cached.cache_id != cur_cache_id_) // = cached version is outdated
      {
        rbi_cached.set_particle(spheres_table,
                                quaternions_tables,
                                pi,
                                cur_cache_id_);
      }
    return rbi_cached;
  }
  else // if is_cache_active_
    {
      return internal::RigidBodyInfo(spheres_table,
                                     quaternions_tables,
                                     pi, INVALID_CACHE_ID);
    }
}



IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_SITES_PAIR_SCORE_H */
