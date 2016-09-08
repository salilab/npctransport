/**
 *  \file IMP/npctransport/SlabWithToroidalPoreSingletonScore.h
 *  \brief a score for a slab with a toroidal pore
 *

 *  Copyright 2007-2015 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPNPCTRANSPORT_SLAB_WITH_TOROIDAL_PORE_SINGLETON_SCORE_H
#define IMPNPCTRANSPORT_SLAB_WITH_TOROIDAL_PORE_SINGLETON_SCORE_H

#include "npctransport_config.h"
#include <IMP/Model.h>
#include <IMP/Pointer.h>
#include <IMP/check_macros.h>
#include <IMP/SingletonScore.h>
#include <IMP/singleton_macros.h>
#include "IMP/core/XYZR.h"

IMPNPCTRANSPORT_BEGIN_NAMESPACE

class IMPNPCTRANSPORTEXPORT SlabWithToroidalPoreSingletonScore : public SingletonScore
{
  double bottom_;  // bottom of slab on x-axis
  double top_;  // top of slab on z-axis
  double midZ_;  // (top + bottom) / 2, for caching some calculations
  double R_;
  double r_;
  double k_;

public:
  //! Constructs a horizontal slab with a toroidal pore,
  //! centered at the z=0 plane
  /**
     Constructs a horizontal slab with a toroidal pore, centered at
     the z=0 plane, from z=-0.5*thickness to z=+0.5*thickness.

      @param radius outer radius of the torus
      @param slab_thickness vertical thickness of slab. The
             minor radius of the pore torus equals half
             the slab thickness.
      @param the slab repulsive force constant in kcal/mol/A
  */
  SlabWithToroidalPoreSingletonScore
    (double slab_thickness, double radius, double k);

  //! Constructs a slab from z=bottom to z=top with a toroidal pore
  /**
      Constructs a slab with a toroidal pore, parallel to the plane z=0
      and running from z=bottom to z=top.
      The minor radius of the pore torus equals half the slab
      thickness, or 0.5*(top-bottom)

      @param slab_bottom bottom z coordinate of slab
      @param slab_top top z coordinate of slab
      @param radius outer radius of the torus
      @param k the slab repulsive force constant in kcal/mol/A
  */
  SlabWithToroidalPoreSingletonScore
    (double slab_bottom, double slab_top, double radius, double k);


  /** returns the lowest slab z coordinate */
  double get_bottom_z() { return bottom_; }

  /** returns the highest slab z coordinate */
  double get_top_z() { return top_; }

  //! evaluate score for particle pi in model m. If da is not null,
  //! use it to accumulate derivatives in model.
  virtual double evaluate_index
    (Model *m,
     ParticleIndex pi,
     DerivativeAccumulator *da) const IMP_OVERRIDE;

  virtual ModelObjectsTemp do_get_inputs(Model *m,
                                         const ParticleIndexes &pis) const
      IMP_OVERRIDE;

  //! evaluate score for particles pis[lower_bound..upper_bound] in
  //! model m. If da is not null, use it to accumulate derivatives
  //! in model.
  virtual double evaluate_indexes
    (Model *m,
     const ParticleIndexes &pis,
     DerivativeAccumulator *da,
     unsigned int lower_bound,
     unsigned int upper_bound) const IMP_FINAL;

  /**
     Evaluate score for particles pis[lower_bound..upper_bound] in
     model m. If da is not null, use it to accumulate derivatives
     in model.

     Abort early and return maximal double value if score is larger
     than max (note this assumes all score components are
     non-negative, which is true for slab score)
  */
  double evaluate_if_good_indexes
    ( Model *m,
      const ParticleIndexes &p,
      DerivativeAccumulator *da,
      double max,
      unsigned int lower_bound,
      unsigned int upper_bound) const
  {
    double ret = 0;
    for (unsigned int i = lower_bound; i < upper_bound; ++i) {
      ret += evaluate_if_good_index(m, p[i], da, max - ret);
      if (ret > max) return std::numeric_limits<double>::max();
    }
    return ret;
  }


 private:
  // return shortest distance needed to move the sphere in order to
  // remove an overlap between the sphere and the slab surface.
  // Return 0 if ball does not penetrate slab.
  //
  // @param s the sphere to evaluate
  // @param out_displacement if not null and the returned score is positive, *out_displacement
  //                         is used to store the normalized displacement vector for removing
  //                         the overlap between the sphere and the slab (or 0,0,0 if no overlap)
  inline double get_sphere_penetration
    (algebra::Sphere3D s,
     algebra::Vector3D* out_displacement) const;


  IMP_OBJECT_METHODS(SlabWithToroidalPoreSingletonScore);
};

//
inline double
SlabWithToroidalPoreSingletonScore::evaluate_index
(Model *m,
 const ParticleIndex pi,
 DerivativeAccumulator *da) const
{
  IMP_OBJECT_LOG;
  IMP::core::XYZR d(m, pi);
  if (!d.get_coordinates_are_optimized()) return false;
  algebra::Vector3D displacement;
  double score=get_sphere_penetration(d.get_sphere(),
                               da ? &displacement : nullptr);
  IMP_LOG_TERSE("sphere " << d << " score " << score);
  if(da && score>0.0){
    algebra::Vector3D derivative_vector = -k_*displacement;
    IMP_LOG_TERSE(" derivative vector " << derivative_vector);
    d.add_to_derivatives(derivative_vector, *da);
  }
  IMP_LOG_TERSE(std::endl);
  return score;
}

//
inline double
SlabWithToroidalPoreSingletonScore::evaluate_indexes
    (Model *m,
     const ParticleIndexes &pis,
     DerivativeAccumulator *da,
     unsigned int lower_bound,
     unsigned int upper_bound) const
{
  IMP_LOG_TERSE("SlabWithToroidalPore singleton - evaluate indexes"
          << std::endl);
  double ret = 0;
  algebra::Sphere3D const* spheres_table=
    m->access_spheres_data();
  algebra::Sphere3D* sphere_derivatives_table=
    m->access_sphere_derivatives_data();
  IMP::internal::BoolAttributeTableTraits::Container const&
    is_optimizable_table= m->access_optimizeds_data
    (core::XYZ::get_coordinate_key(0)); // x is indicator
  // Evaluate and sum score and derivative for all particles:
  for (unsigned int i = lower_bound; i < upper_bound; ++i) {
    int pi_index=pis[i].get_index();
    // Check attributes have valid values:
    IMP_CHECK_CODE( {
        IMP::core::XYZR d(m, pis[i]);
        algebra::Sphere3D s=spheres_table[pi_index];
        IMP_INTERNAL_CHECK(d.get_coordinates_are_optimized() == is_optimizable_table[pis[i]],
                           "optimable table inconsistent with d.get_coordinates_are_optimized for particle " << d);
        IMP_INTERNAL_CHECK((d.get_coordinates() - s.get_center()).get_magnitude()<.001,
                           "Different coords for particle " << d << " *** "
                           << d.get_coordinates() << " vs. " << s.get_center());
        IMP_INTERNAL_CHECK(d.get_radius() == s.get_radius(),
                           "Different radii for particle " << d << " *** "
                           << d.get_radius() << " vs. " << s.get_radius());
      } ); // IMP_CHECK_CODE
    if(!is_optimizable_table[pis[i]]) {
      continue;
    }
    algebra::Vector3D displacement;
    double cur_score = get_sphere_penetration(spheres_table[pi_index],
                                        da ? &displacement : nullptr);
    IMP_LOG_TERSE("SlabWithToroidalPore singleton score for sphere / displacement "
            << spheres_table[pi_index] << " is " << cur_score << " / "
            << displacement << std::endl);
    ret+= cur_score;
    if(cur_score>0.0 && da) {
      algebra::Vector3D derivative_vector = -k_*displacement;
      // accumulate derivatives directly for speed
      for(unsigned int j=0; j<3; j++) {
        sphere_derivatives_table[pi_index][j] += (*da)(derivative_vector[j]);
      }
    }
  }

  return ret;
}

// return - distance to nearest surface point
// out_displacement - a vector pointing to the nearest surface point
double
SlabWithToroidalPoreSingletonScore::get_sphere_penetration
(algebra::Sphere3D sphere,
 algebra::Vector3D* out_displacement) const
{
  double const x=sphere[0];
  double const y=sphere[1];
  double const z=sphere[2];
  double const sr=sphere.get_radius();
  double sphere_dz_top= (z - sr) - top_;
  double sphere_dz_bottom = (z + sr) - bottom_;
  bool is_above_top=sphere_dz_top>0;
  bool is_below_bottom=sphere_dz_bottom<0;
  if(is_above_top || is_below_bottom) {
    // early abort I - above or below slab
    if(out_displacement){
      (*out_displacement)=IMP::algebra::Vector3D(0,0,0);
    }
    return 0;
  }
  double d_xy2 = x*x + y*y; // distance from (0,0,z)
  bool is_closer_to_top= (sphere_dz_bottom > -sphere_dz_top);
  double abs_sphere_dz =
    is_closer_to_top ? -sphere_dz_top : sphere_dz_bottom; // d_z is the mianimal vertical change in z that brings s either fully above or fully below the slab
  if (d_xy2 > R_*R_) {
      // early about II - radially outside torus major radius R
    if(out_displacement){
      (*out_displacement)=IMP::algebra::Vector3D(0,0,is_closer_to_top?+1:-1);
    }
    return abs_sphere_dz;
  }
  double d_xy = std::sqrt(d_xy2);
  if(d_xy+sr < R_-r_){
    // early abort III - radially close to central axis (less than R-r)
    if(out_displacement){
      (*out_displacement)=IMP::algebra::Vector3D(0,0,0);
    }
    return 0;
  }
  // IV. In slab vertically + radially between R-r and R relative to central axis
  const double eps = 1e-9;
  double dR_x, dR_y;
  if ( d_xy > eps )
  {
    // (dR_x,dR_y) is the vector from the nearest point on the center line of the torus to (x,y) (projected on z=0)
    dR_x = x - x*(R_/d_xy);
    dR_y = y - y*(R_/d_xy);
  }
  else
  {
    dR_x = x - R_;
    dR_y = y;
  }
  double dR_z = z-midZ_;
  double dR = std::sqrt(dR_z*dR_z + dR_x*dR_x + dR_y*dR_y); // magnitude of vector from nearest point on the torus central line (radius R circle on z=0)
  double sphere_penetration= r_ + sphere.get_radius() - dR;
  // No overlap:
  if(sphere_penetration<=0.0){
    if(out_displacement){
      (*out_displacement)=IMP::algebra::Vector3D(0,0,0);
    }
    return 0;
  }
  // Overlap:
  if(out_displacement){
    if ( dR > eps )
      {
        (*out_displacement)= IMP::algebra::Vector3D(dR_x,dR_y,dR_z)/dR;
      }
    else
      {
        (*out_displacement)= IMP::algebra::Vector3D(dR_x,dR_y,dR_z)/eps;
      }
  }
  return sphere_penetration;
}



IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_SLAB_WITH_TOROIDAL_PORE_SINGLETON_SCORE_H */
