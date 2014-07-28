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
#include "SitesPairScoreParams.h"

#include <IMP/core/rigid_bodies.h>
#include <IMP/base/log_macros.h>
#include <IMP/algebra/vector_generators.h>

#include <cmath> // overloaded versions of abs!!

IMPNPCTRANSPORT_BEGIN_INTERNAL_NAMESPACE


inline double evaluate_one_site(
    double k, double range, core::RigidBody rb0, core::RigidBody rb1,
    algebra::Transformation3D tr0, algebra::Transformation3D tr1,
    algebra::Rotation3D irot0, algebra::Rotation3D irot1,
    const algebra::Vector3D &local0, const algebra::Vector3D &local1,
    DerivativeAccumulator *da) {
  algebra::Vector3D g0 = tr0.get_transformed(local0);
  algebra::Vector3D g1 = tr1.get_transformed(local1);
  algebra::VectorD<3> delta = g0 - g1;
  static const double MIN_DISTANCE = .00001;
  double distance2 = delta.get_squared_magnitude();
  if (distance2 > square(range)) return 0;
  double distance = std::sqrt(distance2);
  double score;
  score = .5 * k * distance2 - .5 * k * square(range);
  if (da && distance > MIN_DISTANCE) {
    algebra::Vector3D duv = k * delta;
    rb0.add_to_derivatives(irot0.get_rotated(duv), duv, local0,
                           tr0.get_rotation(), *da);
    rb1.add_to_derivatives(irot1.get_rotated(-duv), -duv, local1,
                           tr1.get_rotation(), *da);
  }
  return score;
}
/*inline float Q_rsqrt( float number )
{
  int i;
  float x2, y;
  const float threehalfs = 1.5F;

  x2 = number * 0.5F;
  y  = number;
  i  = *reinterpret_cast<int*>(&y);
  i  = 0x5f3759df - ( i >> 1 );
  y  = *reinterpret_cast<float*>(&i);
  y  = y * ( threehalfs - ( x2 * y * y ) );
  y  = y * ( threehalfs - ( x2 * y * y ) );

  return y;
  }*/

/**
   Evaluates the attraction score between two sites (local0 and local1),
   transformed by tr0 and tr1 respectively, based on their distance

   The returned score is linear -k*(2-distance) for distance < 1.0
   or -k/distance for distance > 1.0
   TODO: why arbitrary threshold distance for linear range set to 1.0?

   @param rb0[out],rb1[out] the particles that contain the sites
   @param tr0, tr1 the particles frames of reference
   @param irot0, irot1 the reverse rotations from each particle frame of ref.
   @param local0,local1 the sites
   @param da[out] derivative accumulator that gets updated
                  with the score derivative

   @return the score
*/
inline double evaluate_one_site_2(
    double k, double range, core::RigidBody rb0, core::RigidBody rb1,
    algebra::Transformation3D tr0, algebra::Transformation3D tr1,
    algebra::Rotation3D irot0, algebra::Rotation3D irot1,
    const algebra::Vector3D &local0, const algebra::Vector3D &local1,
    DerivativeAccumulator *da) {
  algebra::Vector3D g0 = tr0.get_transformed(local0);
  algebra::Vector3D g1 = tr1.get_transformed(local1);
  algebra::VectorD<3> delta = g0 - g1;
  static const double MIN_DISTANCE = .00001;
  double distance2 = delta.get_squared_magnitude();
  if (distance2 > algebra::get_squared(range)) return 0;
  if (distance2 < MIN_DISTANCE) return -k * range;
  // double distance=std::sqrt(distance2);
  // double dp1= distance+1;
  // double dp12= algebra::get_squar4ed(dp1);
  double idistance = 1.0f / sqrtf(distance2);
  double distance = distance2 * idistance;
  double kidistance = k * idistance;
  algebra::VectorD<3> deriv = kidistance * delta;  // magnitude k
  double score = -k * (range - distance);
  if (da && distance2 > MIN_DISTANCE) {
    rb0.add_to_derivatives(irot0.get_rotated(deriv), deriv, local0,
                           tr0.get_rotation(), *da);
    rb1.add_to_derivatives(irot1.get_rotated(-deriv), -deriv, local1,
                           tr1.get_rotation(), *da);
  }
  return score;
}

// note it can be negative if spheres overlap so not strictly
// 'distance' (in this case the negative result is the distance needed
// to make the spheres touch)
//
// gD = gD*[1-(r1+r2)/|gD|]
//
// @param gD [in/out] - raw distance vector to be corrected
// @param r0 - radius of sphere 0
// @param r1 - radius of sphere 1
//
// @return the magnitude of the corrected gD
inline double subtract_sphere_radii_from_distance_vector
(algebra::Vector3D& gD, double r0, double r1)
{
  static double const MIN_D = 0.001;
  double d2 = gD.get_squared_magnitude();
  double d = sqrtf(d2);
  double dr= d - r0 - r1;
  if(d > MIN_D){
    gD *= (dr/d);
    return dr;
  } else {
    gD = dr * IMP::algebra::get_random_vector_on_unit_sphere();
  }
  return dr;
}


/**
   Evaluate interaction on a pair of sites as a linear potential
   (constant force) within the passed attraction range.

   The maximal drop in potential energy in kcal/mol is:
     DELTA-U = k*range

    @param k - force constant
    @param range - attraction range
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


// dX - distance (assumed non-negative)
// spsp
// f - force magnitude (output; always non-negative if k positive)
//
// returns energy potential (always non-positive if k positive)
inline double get_U(double dX,
                         SitesPairScoreParams const& spsp,
                         double& f)
{
  double dX2 = dX*dX;
  double const& r = spsp.r;
  double const& k = spsp.k;
  double const& kr = spsp.kr;
  double const& kr2 = spsp.kr2;
  //  int sign = dX > 0 ? 1 : -1;
  //dX *= sign; // make non-negative
  if(dX<r/2)
    {
      f= k*dX;//*sign;
      return 0.5*k*dX2-0.25*kr2;
    }
  if(dX<r)
    {
      f=k*(r-dX);//*sign;
      return -0.5*k*dX2 + kr*dX - 0.5*kr2;
    }
  f = 0;
  return 0;
}

// dX,dY - distance in x,y directions (assumed non-negative!)
// params_X, paramsY
// fX,fY - force magnitude in x,y directions (output; always non-negative if k positive)
//
// returns energy potential (always non-positive if k positives)
inline double get_V(double dX, double dY,
                         SitesPairScoreParams const& params_X,
                         SitesPairScoreParams const& params_Y,
                         double &fX, double& fY)
{
  double fX1, fY1;
  IMP_LOG_PROGRESS(dX << ","
                   << dY << ","
                   << params_X.k << ","
                   << params_X.r << ","
                   << params_Y.k << ","
                   << params_Y.r << ","
                   << std::endl);
  double UX=get_U(dX, params_X, fX1);
  double UY=get_U(dY, params_Y, fY1);
  fX = -UY * fX1;
  fY = -UX * fY1;
  return -UX * UY;
}


//! Evaluate interaction on pair of sites
/** Evaluate interaction on pair of sits using a bell-shaped spline
    function with force coefficient k and range r, which is skewed
    (anisotropic) in normal vs. tangent direction, with normal
    direction defined as the vector between the two beads of rbi0 and
    rbi1/

    @param params_unskewed params bounds regardless of directionality
    @param params_N range and k params in normal direction
    @param params_T range and k params in tangent direction
    @param rbi0 - cached information on rigid body 0
    @param rbi1 - cached information on rigid body 1
    @param l0 - site0 local coordinates + radius
    @param l1 - site1 local coordinates + radius
    @param g0 - site0 global coordinates
    @param g1 - site1 global coordinates
    @param da - accumulator for reweighting derivatives
    */
inline double evaluate_one_site_4
( SitesPairScoreParams const& params_unskewed,
  SitesPairScoreParams const& params_N,
  SitesPairScoreParams const& params_T,
  RigidBodyInfo& rbi0, RigidBodyInfo& rbi1,
  const algebra::Sphere3D &l0, const algebra::Sphere3D &l1,
  algebra::Vector3D& g0, algebra::Vector3D& g1,
  DerivativeAccumulator *da)
{
  double const& range = params_unskewed.r;
  algebra::Vector3D gD = g0 - g1; // g = global
  // Quick filters:
  static const double MIN_D = 0.001;
  IMP_LOG_PROGRESS("point distance = " << sqrtf((g0-g1).get_squared_magnitude()));
  double d = subtract_sphere_radii_from_distance_vector(gD, l0.get_radius(), l1.get_radius());
  IMP_LOG_PROGRESS(" ; r0= " << l0.get_radius()
                   << ", ; r1 = " << l1.get_radius()
                   << " ; sphere distance = " << d
                   << " ; MIN_D = " << MIN_D
                   << " ; range = " << range
                   << std::endl);
  if (d > range){ // (d2 > range2)
    IMP_LOG_PROGRESS( "beyond range" << std::endl);
    return 0;
  }
  if (std::abs(d) < MIN_D) {
    IMP_LOG_PROGRESS( "d < MIN_D" << std::endl);
    //    return -0.25 * params_unskewed.k * params_unskewed.r2;
    return 0.0625 * params_unskewed.k * params_N.r2 * params_T.r2;
  }
  // Deconvolute gD to normal and tangent components
  // gD = dN * gN + dT * gT
  // where normal direction = vector connecting two rigid bodies:
  double d2 = d*d;
  algebra::VectorD<3> gN =
    rbi0.rb.get_coordinates() - rbi1.rb.get_coordinates();
  gN /= gN.get_magnitude();
  double dN = gN * gD; // note: may be negative
  double dT = sqrtf(d2 - std::pow(dN,2)); // note: non-negative by constr.
  double score, Fx, Fy;
  score = get_V(std::abs(dN), dT, params_N, params_T, Fx, Fy);
  IMP_LOG_PROGRESS("site score " << score << std::endl);
  IMP_LOG_PROGRESS( " gD " << gD
                    << " ; gN " << gN
                    << " ; dN " << dN
                    << " ; dT " << dT
                    << std::endl);

  if (da)
    {
      algebra::Vector3D gderiv0(0,0,0);
      if(Fx > 0) {
        gderiv0 = Fx * gN;
        if(dN < 0){ gderiv0 *= -1.0; }
      }
      if(Fy > 0) {
        algebra::Vector3D gT = (gD - dN * gN) / dT;
        gderiv0 += Fy * gT;
      IMP_LOG(PROGRESS,  "gT " << gT);
      }
      IMP_LOG(PROGRESS, "Fx " << Fx << " ; Fy " << Fy
              << " ; gderiv0 " << gderiv0
              << std::endl);
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


IMPNPCTRANSPORT_END_INTERNAL_NAMESPACE

#endif /* IMPNPCTRANSPORT_INTERNAL_SITES_H */
