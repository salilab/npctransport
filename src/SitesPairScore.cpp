/**
 *  \file DistancePairScore.cpp
 *  \brief A Score on the distance between a pair of particles.
 *
 *  Copyright 2007-2018 IMP Inventors. All rights reserved.
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
  sites1_(sites1) // ,
  //  is_cache_active_(false),
  //  cur_cache_id_(INVALID_CACHE_ID)
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

//!
double
SitesPairScore::evaluate_indexes
(Model *m, const ParticleIndexPairs &pis,
 DerivativeAccumulator *da,
 unsigned int lower_bound,
 unsigned int upper_bound) const
{
  // get internal tables:
  algebra::Sphere3D const* spheres_table=
      m->access_spheres_data();
  algebra::Sphere3D* sphere_derivatives_table=
    m->access_sphere_derivatives_data();
  double const* quaternions_tables[4];
  for(unsigned int i = 0; i < 4; i++){
    quaternions_tables[i]=
      core::RigidBody::access_quaternion_i_data(m, i);
  }
  double* torques_tables[3];
  for(unsigned int i = 0; i < 3; i++){
    torques_tables[i]=
      core::RigidBody::access_torque_i_data(m, i);
  }
  // evaluate all idexes with rigid body info cache active:
  //   activate_cache();
  double ret = 0.0;
  for (unsigned int i = lower_bound; i < upper_bound; ++i) {
    ret += evaluate_index_with_internal_tables(m,
                                               spheres_table,
                                               quaternions_tables,
                                               sphere_derivatives_table,
                                               torques_tables,
                                               pis[i],
                                               da);
  }
  //   deactivate_cache();
  return ret;
}



ModelObjectsTemp
SitesPairScore::do_get_inputs(
    Model *m, const ParticleIndexes &pis) const {
  return IMP::get_particles(m, pis);
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
