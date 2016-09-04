/**
 *  \file IMP/npctransport/SlabWithToroidalPoreSingletonScore.h
 *  \brief a score for a slab with a toroidal pore
 *

 *  Copyright 2007-2015 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPNPCTR
#define IMPNPCTRANSPORT_SLAB_WITH_TORUS_SINGLETON_SCORE_H

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
      @param thickness vertical thickness of slab. The
             minor radius of the pore torus equals half
             the slab thickness.
      @param the slab repulsive force constant in kcal/mol/A
  */
  SlabWithToroidalPoreSingletonScore
    (double radius, double thickness, double k);

  //! Constructs a slab from z=bottom to z=top with a toroidal pore
  /**
      Constructs a slab with a toroidal pore, parallel to the plane z=0
      and running from z=bottom to z=top.
      The minor radius of the pore torus equals half the slab
      thickness, or 0.5*(top-bottom)

      @param radius outer radius of the torus
      @param bottom bottom z coordinate of slab
      @param top top z coordinate of slab
      @param the slab repulsive force constant in kcal/mol/A
  */
  SlabWithToroidalPoreSingletonScore
    (double radius, double bottom, double top, double k);


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
  // evaluate slab for specified sphere. Return 0 if ball
  // does not penetrate slab
  //
  // @param s the sphere to evaluate
  // @param out_displacement if not null and the returned score is positive, *out_displacement
  //                         is used to store the computed displacement vector from
  //                         the surface of the z-axis aligned cylinder to
  //                         the center of s. Ignore if score is zero.
  inline double evaluate_sphere
    (algebra::Sphere3D s,
     algebra::Vector3D* out_displacement) const;


  // computes the displacement from the surface of the z-axis
  // aligned cylinder to v
  //
  // @return <distance, a vector pointing out>,
  //         negative distance means v is inside cylinder
  std::pair<double, algebra::Vector3D> get_displacement_vector(
      const algebra::Vector3D &v) const;

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
  double score=evaluate_sphere(d.get_sphere(),
                               da ? &displacement : nullptr);
  if(da && score>0.0){
    algebra::Vector3D derivative_vector = -k_*displacement;
    IMP_LOG(PROGRESS, "result in " << score << " and " << derivative_vector << std::endl);
    d.add_to_derivatives(derivative_vector, *da);
  }
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
    // Check attributes have valid valies:
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
    double cur_score = evaluate_sphere(spheres_table[pi_index],
                                        da ? &displacement : nullptr);
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

//
double SlabWithToroidalPoreSingletonScore::evaluate_sphere
(algebra::Sphere3D sphere,
 algebra::Vector3D* out_displacement) const
{
  double const x=sphere[0];
  double const y=sphere[1];
  double const z=sphere[2];
  double const sr=sphere.get_radius();
  double dz_top= (z - sr) - top_;
  double dz_bottom = bottom_ - (z + sr);
  if ((dz_top > 0) || (dz_bottom>0)) {
    // early abort I - above or below slab
    if(out_displacement){
      (*out_displacement)=IMP::algebra::Vector3D(0,0,0);
    }
    return 0;
  }
  double d_xy2 = x*x + y*y; // distance from (0,0,z)
  double d_z = (dz_top < dz_bottom) ? -dz_top : dz_bottom; // d_z is the minimal vertical change in z that brings s either fully above or fully below the slab
  if (d_xy2 > R_*R_) {
      // early about II - outside torus center line
    if(out_displacement){
      (*out_displacement)=IMP::algebra::Vector3D(0,0,(d_z>0)?1:-1);
    }
    return d_z;
  }
  double d_xy = std::sqrt(d_xy2);
  if(d_xy+sr < R_-r_){
    // early abort III - deep within torus hole
    if(out_displacement){
      (*out_displacement)=IMP::algebra::Vector3D(0,0,0);
    }
    return 0;
  }
  const double eps = 1e-9;
  double d_tx, d_ty;
  if ( d_xy > eps )
  {
    // (d_tx,d_ty) is the vector from the nearest point on the center line of the torus to (x,y) (projected on z=0)
    d_tx = x - x/d_xy*R_;
    d_ty = y - y/d_xy*R_;
  }
  else
  {
    d_tx = x - R_;
    d_ty = y;
  }
  double d_tz = z-midZ_;
  double denom = std::sqrt(d_tz*d_tz + d_tx*d_tx + d_ty*d_ty); // magnitude of vector from nearest point on the torus center line to
  if(out_displacement){
    if ( denom > eps )
      {
        (*out_displacement)= IMP::algebra::Vector3D(d_tx,d_ty,d_tz)/denom;
      }
    else
      {
        (*out_displacement)= IMP::algebra::Vector3D(d_tx,d_ty,d_tz)/eps;
      }
  }
  return denom - r_ - sphere.get_radius();
}



IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_SLAB_WITH_TORUS_SINGLETON_SCORE_H */
