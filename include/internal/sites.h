/**
 *  \file SitesPairScore.h
 *  \brief A Score on the distance between a pair of particles.
 *
 *  Copyright 2007-2010 IMP Inventors. All rights reserved.
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
  DerivativeAccumulator *da)
{
  static const double MIN_D = .001;
  algebra::VectorD<3> gD = g0 - g1;
  double d = subtract_sphere_radii_from_distance_vector(gD, l0.get_radius(), l1.get_radius());
  d = std::abs(d);
  if (d > range) return 0;
  if (d < MIN_D) return -k * range;
  double id = 1.0f / d;
  //  double d = d2 * id;
  double score = -k * (range - d);
  if (da) {
    algebra::Vector3D gderiv0 = (k * id) * gD;
    algebra::Vector3D lderiv0 = rbi0.irot.get_rotated(gderiv0);
    rbi0.rb.add_to_derivatives(lderiv0, gderiv0, l0.get_center(),
                               rbi0.tr.get_rotation(),*da);
    algebra::Vector3D gderiv1 = -gderiv0;
    algebra::Vector3D lderiv1 = rbi1.irot.get_rotated(gderiv1);
    rbi1.rb.add_to_derivatives(lderiv1, gderiv1, l1.get_center(),
                               rbi1.tr.get_rotation(),*da);
  }
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
    @param lSite1 - site0 local coordinates + radius
    @param lSite2 - site1 local coordinates + radius
    @param gSite1 - site0 global coordinates
    @param gSite2 - site1 global coordinates
    @param da - accumulator for reweighting derivatives,
                or null to disable force and torque computations
    */
inline
double evaluate_pair_of_sites
( SitesPairScoreParameters const& spsp,
  RigidBodyInfo& rbi1, RigidBodyInfo& rbi2,
  const algebra::Sphere3D &lSite1, const algebra::Sphere3D &lSite2,
  algebra::Vector3D& gSite1, algebra::Vector3D& gSite2,
  DerivativeAccumulator *da)
{
  using IMP::algebra::Vector3D;
  using IMP::algebra::Sphere3D;

  // I. Pre-computations:
  Vector3D gRB1=rbi1.rb.get_coordinates();
  Vector3D gRB2=rbi2.rb.get_coordinates();
 // update coordinates
  Vector3D gRB1RB2 = gRB2-gRB1;
  Vector3D gUnitRB1RB2 = gRB1RB2.get_unit_vector();
  Vector3D gUnitRB2RB1 = -gUnitRB1RB2;
  Vector3D gUnitRB1Site1 = (gSite1-gRB1)*rbi1.iradius;
  Vector3D gUnitRB2Site2 = (gSite2-gRB2)*rbi2.iradius;
  double cosSigma1 = gUnitRB1Site1*gUnitRB1RB2;
  double cosSigma2 = gUnitRB2Site2*gUnitRB2RB1;
  if(cosSigma1<spsp.cosSigma1_max || cosSigma2<spsp.cosSigma2_max){
    // equivalent to sigma>sigma_max, so out of range
    return 0;
  }

  // II. Energy computations:
  double kFactor1=get_k_factor(cosSigma1, spsp.cosSigma1_max);
  double kFactor2=get_k_factor(cosSigma2, spsp.cosSigma2_max);
  double kFactor=kFactor1*kFactor2;
  IMP_LOG_VERBOSE('kFactor1 ' << kFactor1 <<
		  ' kFactor2 ' << kFactor2);
  double u_1D, derivR_1D; // energy and derivative before factoring k
  double rbSphereDist=gRB1RB2.get_magnitude()-rbi1.radius-rbi2.radius;
  u_1D=get_U_1D(rbSphereDist, spsp, derivR_1D);
  double score=kFactor*u_1D;
  IMP_LOG_VERBOSE('score ' << score);

  // III. Apply force to rigid bodies if derivative accumulator is active
  if(da){
    // Add translational force:
    double derivR=kFactor*derivR_1D;
    Vector3D gDerivR_on_RB1=derivR*gUnitRB1RB2;
    IMP::core::XYZ(rbi1.rb).add_to_derivatives( gDerivR_on_RB1, *da);  // action
    IMP::core::XYZ(rbi2.rb).add_to_derivatives(-gDerivR_on_RB1, *da); // reaction
    IMP_LOG_VERBOSE('global translation derivative on first rb ' << gDerivR_on_RB1);
    // Add torque:
    // (note it is assumed that the opposing torque is
    // dissipated in water, so no action/reaction between RB1 and RB2)
    if(kFactor1>0.0 && kFactor1<0.99999){ // within attraction range
      Vector3D tmp1 = get_vector_product(gUnitRB1Site1,gUnitRB1RB2);
      double absSinSigma1 = tmp1.get_magnitude();
      IMP_USAGE_CHECK(absSinSigma1>0.00001,
		      "abs(sinSigma) is expected to be positive within attraction force range" );
      Vector3D gRotSigma1=tmp1/absSinSigma1; // rotation axis about RB1 that brings site1 towards the RB1-RB2 axis
      double dKFactor1=get_derivative_k_factor(absSinSigma1, spsp.cosSigma1_max);
      double fS1=-u_1D*dKFactor1*kFactor2;
      Vector3D gTorque_on_RB1=fS1*gRotSigma1;
      Vector3D lTorque_on_RB1=rbi1.irot.get_rotated(gTorque_on_RB1);
      rbi1.rb.add_to_torque(lTorque_on_RB1, *da);
      IMP_LOG_VERBOSE('global torque on first rb ' << lTorque_on_RB1);
    }
    if(kFactor2>0.0 && kFactor2<0.99999){ // within attraction range
      Vector3D tmp2 = get_vector_product(gUnitRB2Site2,gUnitRB2RB1);
      double absSinSigma2 = tmp2.get_magnitude();
      IMP_USAGE_CHECK(absSinSigma2>0.00001,
		      "abs(sinSigma2) is expected to be positive within attraction force range" );
      Vector3D gRotSigma2=tmp2/absSinSigma2; // rotation axis about RB2 that brings site2 towards the RB2-RB1 axis
      double dKFactor2=get_derivative_k_factor(absSinSigma2, spsp.cosSigma2_max);
      double fS2=-u_1D*kFactor1 *dKFactor2;
      Vector3D gTorque_on_RB2=fS2*gRotSigma2;
      Vector3D lTorque_on_RB2=rbi2.irot.get_rotated(gTorque_on_RB2);
      rbi2.rb.add_to_torque(lTorque_on_RB2, *da);
      IMP_LOG_VERBOSE('global torque on second rb ' << lTorque_on_RB2);
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
