/**
 *  \file SitesPairScore.h
 *  \brief A Score on the distance between a pair of particles.
 *
 *  Copyright 2007-2017 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_INTERNAL_SITES_H
#define IMPNPCTRANSPORT_INTERNAL_SITES_H

#include "../npctransport_config.h"
#include "RigidBodyInfo.h"
#include "../SitesPairScoreParameters.h"

#include <IMP/core/rigid_bodies.h>
#include <IMP/log_macros.h>
#include <IMP/algebra/vector_generators.h>

#include <cmath> // overloaded versions of abs!!

IMPNPCTRANSPORT_BEGIN_INTERNAL_NAMESPACE

// return the sphere distance between the end points of gD,
// and overrides gD with the distance vector between the corresponding
// sphere surfaces
//
// gD = gD*[1-(r1+r2)/|gD|]
//
// @param gD [in/out] - distance vector between sphere centers,
//                      to be corrected to the distance vector
//                      between the sphere surfaces
// @param r0 - radius of sphere 0
// @param r1 - radius of sphere 1
//
// @return the magnitude of the corrected gD
inline double subtract_sphere_radii_from_distance_vector
(algebra::Vector3D& gD, double r0, double r1)
{
  static double const MIN_D = 0.001;
  double d = gD.get_magnitude();
  double dr= d - r0 - r1;

  if(d > MIN_D){
    gD *= (dr/d);
  } else {
    gD = dr * IMP::algebra::get_random_vector_on_unit_sphere();
  }
  return std::abs(dr);
}


/**
   Evaluate interaction on a pair of sites as a linear potential
   (constant force) within the passed attraction range.

   The maximal drop in potential energy in kcal/mol is:
     DELTA-U = k*range

    @param k - force constant [kCal/mol/A]
    @param range - attraction range [A]
    @param rbi0 - cached information on rigid body 0
    @param rbi1 - cached information on rigid body 1
    @param l0 - site0 local coordinates + radius
    @param l1 - site1 local coordinates + radius
    @param g0 - site0 global coordinates
    @param g1 - site1 global coordinates
    @param da - accumulator for reweighting derivatives
    */
inline double evaluate_one_site_3
( double k,
  double range,
  RigidBodyInfo& rbi0, RigidBodyInfo& rbi1,
  const algebra::Sphere3D &l0, const algebra::Sphere3D &l1,
  algebra::Vector3D& g0, algebra::Vector3D& g1,
  DerivativeAccumulator *da,
  algebra::Sphere3D *sphere_derivatives_table,
  double **torques_tables)
{
  static const double MIN_D = .001;
  algebra::Vector3D gD = g0 - g1;
  double d = subtract_sphere_radii_from_distance_vector(gD, l0.get_radius(), l1.get_radius());
  d = std::abs(d);
  if (d > range) return 0;
  if (d < MIN_D) return -k * range;
  double id = 1.0f / d;
  //  double d = d2 * id;
  double score = -k * (range - d);
  if (da) {
    // RB0:
    algebra::Vector3D gderiv0 = (k * id) * gD;
    algebra::Vector3D lderiv0 = rbi0.irot.get_rotated(gderiv0);
    //core::RigidBody(rbi0.pi).add_to_derivatives(lderiv0, gderiv0, l0.get_center(),
    //                                           rbi0.tr.get_rotation(),*da); // TODO: swith to internal tables access?
    algebra::Vector3D torque0 = algebra::get_vector_product(l0.get_center(), lderiv0);
    for (unsigned int i = 0; i < 3; ++i) {
      sphere_derivatives_table[rbi0.pi.get_index()][i]+= (*da)(gderiv0[i]);
      torques_tables[i][rbi0.pi.get_index()]+= (*da)(torque0[i]);
    }
    // RB1:
    algebra::Vector3D gderiv1 = -gderiv0;
    algebra::Vector3D lderiv1 = rbi1.irot.get_rotated(gderiv1);
    //    core::RigidBody(rbi1.pi).add_to_derivatives(lderiv1, gderiv1, l1.get_center(),
    //                                           rbi1.tr.get_rotation(),*da);
    algebra::Vector3D torque1 = algebra::get_vector_product(l1.get_center(), lderiv1);
    for (unsigned int i = 0; i < 3; ++i) {
      sphere_derivatives_table[rbi1.pi.get_index()][i]+= (*da)(gderiv1[i]);
      torques_tables[i][rbi1.pi.get_index()]+= (*da)(torque1[i]);
    }
  } // if (da)
  return score;
}




/**
  Computes a bell-shaped attractive potential with range spsp.r,
  which results from an attractive force that increases linearly
  from dX=0 to dX=0.5*spsp.r with slope spsp.k, and diminishes linearly
  between dX=0.5*spsp.r and dX=spsp.r with slope -spsp.k.

  @param dX distnace (note negative dX mean overlap, which will result
                      in spsp.k*dX repulsion force)
  @param spsp score parameters and precomputed parameters, including
              range spsp.r in A, and force coefficient spsp.k in units of
              kCal/mol/A.
  @param deriv derivative of U_1D at dX, opposite sign to spsp.k*dX)
*/
inline double get_U_1D(double dX,
                         SitesPairScoreParameters const& spsp,
                         double& deriv)
{
  double dX2 = dX*dX;
  double const& r = spsp.r;
  double const& k = spsp.k;
  double const& kr = spsp.kr; // k*r
  double const& kr2 = spsp.kr2; // k*r^2
  if(dX<r/2)
    {
      deriv= -k*dX;//*sign;
      return 0.5*k*dX2-0.25*kr2;
    }
  if(dX<r)
    {
      deriv=-k*(r-dX);//*sign;
      return -0.5*k*dX2 + kr*dX - 0.5*kr2; // integral of f, shifted such that U=0 at dX=r.
    }
  // otherwise:
  deriv = 0;
  return 0;
}

inline double get_k_factor(double cos_sigma, double cos_sigma_max){
  return (cos_sigma-cos_sigma_max)/(1.0-cos_sigma_max);
}

inline double get_derivative_k_factor(double sin_sigma, double cos_sigma_max){
  return sin_sigma/(cos_sigma_max-1.0);
}

//! Evaluate anisotropic interaction potential between pair of sites
/** Evaluate interaction score on pair of sites that results from an
    anisotropic attractive translational force, and optionally,
    accumulate appropriate translational and torque derivatives in da.

    The attractive translation force acts on the axis between the
    centers of masses of rbi1.rb and rbi2.rb, it increases linearly
    till mid-range and diminishes linearly between mid-range and
    full-range; see get_U_1D(). The force decays linearly with the
    angles sigma1 and sigma2 formed between the axis of translation
    and either site1 or site2 within a range sigma1_max or sigma2_max,
    respectively.  .  The resulting energy potential results in
    appropriate torques on rb1 and rb2.


bell-shaped spline
    function with force coefficient k and range r, which is skewed
    (anisotropic) in normal vs. tangent direction, with normal
    direction defined as the vector between the two beads of rbi0 and
    rbi1/

    @param spsp parameters of sites pair score
    @param rbi1 - cached information on rigid body 1
    @param rbi2 - cached information on rigid body 2
    @param gSite2 - site1 global coordinates
    @param gUnitRB1RB2 - precomputed unit vector pointing from RB1 towards RB2, in global coordinates
    @param distRB1RB2 - precomputed distance between centers of rb1 and rb2
    @param gRotSigma1 -  precomputed rotation axis of sigma1
    @param kFactor1 - precomputed scaling factor due to sigma1 vs. spsp.cosSigma1_max
    @param dKFactor1 - the pre-computed sigma1 derivative of kFactor1
    @param da - accumulator for reweighting derivatives,
                or null to disable force and torque computations
    @param sphere_derivatives_table
    @param torques_tables
    */
inline
double evaluate_pair_of_sites
( SitesPairScoreParameters const& spsp,
  RigidBodyInfo const& rbi1, RigidBodyInfo const& rbi2,
  algebra::Vector3D const& gSite2,
  algebra::Vector3D const& gUnitRB1RB2,
  double distRB1RB2,
  algebra::Vector3D const& gRotSigma1,
  double kFactor1, double dKFactor1,
  DerivativeAccumulator *da,
  algebra::Sphere3D *sphere_derivatives_table,
  double **torques_tables)
{
  using IMP::algebra::Vector3D;
  using IMP::algebra::Sphere3D;

  // I. Pre-computations:
  Vector3D const& gRB2= rbi2.tr.get_translation();
  Vector3D gUnitRB2RB1= -gUnitRB1RB2;
  Vector3D gUnitRB2Site2 = (gSite2-gRB2)*rbi2.iradius;
  double cosSigma2 = gUnitRB2Site2*gUnitRB2RB1;
  if(cosSigma2<spsp.cosSigma2_max){
    // equivalent to sigma>sigma_max, so out of range
    return 0;
  }

  // II. Energy computations:
  double kFactor2=get_k_factor(cosSigma2, spsp.cosSigma2_max);
  double kFactor=kFactor1*kFactor2;
  IMP_LOG_VERBOSE("kFactor1 " << kFactor1 <<
		  " kFactor2 " << kFactor2);
  double u_1D, derivR_1D; // energy and derivative before factoring k
  double rbSphereDist=distRB1RB2-rbi1.radius-rbi2.radius;
  u_1D=get_U_1D(rbSphereDist, spsp, derivR_1D);
  double score=kFactor*u_1D;
  IMP_LOG_VERBOSE("score " << score);

  // III. Apply force to rigid bodies if derivative accumulator is active
  if(da){
    // Add translational force:
    double derivR=kFactor*derivR_1D;
    Vector3D gDerivR_on_RB1=derivR*gUnitRB1RB2;
    //    IMP::core::XYZ(rbi1.rb).add_to_derivatives( gDerivR_on_RB1, *da);  // action
    //bin/ IMP::core::XYZ(rbi2.rb).add_to_derivatives(-gDerivR_on_RB1, *da); // reaction
    for(unsigned int i=0; i<3; i++){
      double gDerivR_on_RB1_i=(*da)(gDerivR_on_RB1[i]);
      sphere_derivatives_table[rbi1.pi.get_index()][i]+= gDerivR_on_RB1_i; // action (+)
      sphere_derivatives_table[rbi2.pi.get_index()][i]-= gDerivR_on_RB1_i; // opposite reaction (0)
    }
    IMP_LOG_VERBOSE("global translation derivative on first rb " << gDerivR_on_RB1);
    // Add torque:
    // (note it is assumed that the opposing torque is
    // dissipated in water, so no action/reaction between RB1 and RB2)
    if(kFactor1>0.0 && kFactor1<0.99999){ // within attraction range
      // Vector3D gRotSigma1 = get_vector_product(gUnitRB1Site1,gUnitRB1RB2);
      // double absSinSigma1 = get_magnitude_and_normalize_in_place(gRotSigma1);
      // double dKFactor1=get_derivative_k_factor(absSinSigma1, spsp.cosSigma1_max);
      double fS1=-u_1D*dKFactor1*kFactor2;
      Vector3D gTorque_on_RB1=fS1*gRotSigma1;
      Vector3D lTorque_on_RB1=rbi1.irot.get_rotated(gTorque_on_RB1);
      //rbi1.rb.add_to_torque(lTorque_on_RB1, *da);
      for(unsigned int i=0; i<3; i++){
        torques_tables[i][rbi1.pi.get_index()]+= (*da)(lTorque_on_RB1[i]);
      }
      IMP_LOG_VERBOSE("global torque on first rb " << lTorque_on_RB1);
    }
    if(kFactor2>0.0 && kFactor2<0.99999){ // within attraction range
      Vector3D gRotSigma2 = get_vector_product(gUnitRB2Site2,gUnitRB2RB1);
      double absSinSigma2 = get_magnitude_and_normalize_in_place(gRotSigma2);
      double dKFactor2=get_derivative_k_factor(absSinSigma2, spsp.cosSigma2_max);
      double fS2=-u_1D*kFactor1 *dKFactor2;
      Vector3D gTorque_on_RB2=fS2*gRotSigma2;
      Vector3D lTorque_on_RB2=rbi2.irot.get_rotated(gTorque_on_RB2);
      //      rbi2.rb.add_to_torque(lTorque_on_RB2, *da);
      for(unsigned int i=0; i<3; i++){
        torques_tables[i][rbi2.pi.get_index()]+= (*da)(lTorque_on_RB2[i]);
      }
      IMP_LOG_VERBOSE("global torque on second rb " << lTorque_on_RB2);
    }
    /* Vector3D lTorque_on_RB1=rbi1.irot.get_rotated(gTorque_on_RB1-gTorque_on_RB2); // action */
    /* Vector3D lTorque_on_RB2=rbi2.irot.get_rotated(gTorque_on_RB2-gTorque_on_RB1); //reaction */
    /* rbi1.rb.add_to_torque(lTorque_on_RB1, *da); */
    /* rbi2.rb.add_to_torque(lTorque_on_RB2, *da); */
  }

  return score;
}





IMPNPCTRANSPORT_END_INTERNAL_NAMESPACE

#endif /* IMPNPCTRANSPORT_INTERNAL_SITES_H */
