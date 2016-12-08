/**
 *  \file IMP/npctransport/SlabWithToroidalPorePairScore.h
 *  \brief a score for a slab with a toroidal pore
 *

 *  Copyright 2007-2015 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPNPCTRANSPORT_SLAB_WITH_TOROIDAL_PORE_PAIR_SCORE_H
#define IMPNPCTRANSPORT_SLAB_WITH_TOROIDAL_PORE_PAIR_SCORE_H

#include "npctransport_config.h"
#include "SlabWithToroidalPore.h"
#include <IMP/Model.h>
#include <IMP/Pointer.h>
#include <IMP/check_macros.h>
#include <IMP/PairScore.h>
#include <IMP/pair_macros.h>
#include "IMP/core/XYZR.h"

IMPNPCTRANSPORT_BEGIN_NAMESPACE

class IMPNPCTRANSPORTEXPORT SlabWithToroidalPorePairScore : public PairScore
{
  double k_; // repulsion constant in kcal/mol/A

  // cached variables for currently evaluated slab (therefore, mutable):
  mutable double top_; // top z - for caching
  mutable double bottom_; // bottom z - for caching
  mutable double midZ_;  // vertical center of slab
  mutable double R_; // torus major radius
  mutable double rh_; // horizontal minor radius (ellipse horizontal semi-axis)
  mutable double rv_; // vertical minor radius (ellipse vertical semi-axis)
  mutable bool is_pore_radius_optimized_;

public:
  //! Constructs a horizontal slab with a toroidal pore,
  //! centered at the z=0 plane
  /**
     Constructs a score over a horizontal slab with a ring toroidal pore

     @param k the slab repulsive force constant in kcal/mol/A
  */
  SlabWithToroidalPorePairScore
    (double k);

 public:
  //! evaluate score for particle pi in model m. If da is not null,
  //! use it to accumulate derivatives in model.
  virtual double evaluate_index
    (Model *m,
     ParticleIndexPair const& pip,
     DerivativeAccumulator *da) const IMP_OVERRIDE;

  virtual ModelObjectsTemp do_get_inputs(Model *m,
                                         const ParticleIndexes &pis) const
      IMP_OVERRIDE;

  //! evaluate score for particles pis[lower_bound..upper_bound] in
  //! model m. If da is not null, use it to accumulate derivatives
  //! in model.
  virtual double evaluate_indexes
    (Model *m,
     const ParticleIndexPairs &pips,
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
      const ParticleIndexPairs &p,
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

  IMP_OBJECT_METHODS(SlabWithToroidalPorePairScore);

 private:

  /**
     Computes the penetration depth between the specified sphere and an
     axis aligned ellipsoid defined by the torus vertical minor radius and
     horizontal minor radius, centered at the specified origin.
     Optionally stores a normalized translation vector in *out_translation

     @param sphere Sphere for which the distance is computed
     @param origin Ellipsoid origin (a point on the major circle of the torus)
     @param out_translation A pointer to the normalized translation
            vector for removing the sphere out of the ellipsoid is
            stored in *out_translation, if it is not a null pointer.
            (0,0,0) is stored for 0 overlap.

     @return The penetration depth or 0 when the sphere does not overlap with
             the ellipsoid.
  */
  inline double
    get_sphere_ellipsoid_penetration_depth
    (algebra::Sphere3D const& sphere,
     algebra::Vector3D const& origin,
     algebra::Vector3D* out_translation) const;

  // returns the penetration depth D of a sphere relative to the porus slab surface,
  // or 0 if the sphere and the slab do not overlap.
  //
  // @param s the sphere to evaluate
  // @param out_translation if not null, *out_translation
  //                         is used to store the output normalized displacement vector for removing
  //                         the sphere out of the slab (or 0,0,0 if no overlap)
  //                         In other words, D*out_translation be the shortest displacement vector
  //                         needed to remove the sphere out of the slab.
  inline double get_sphere_penetration_depth
    (algebra::Sphere3D s,
     algebra::Vector3D* out_translation) const;

  // update internal variables holding slab params for fast access
  // based on decorated particle swp
  void update_cached_slab_params
    (SlabWithToroidalPore slab) const;
};

inline void
SlabWithToroidalPorePairScore::update_cached_slab_params
(SlabWithToroidalPore slab) const
{
  //TODO: support slabs with non-zero x,y,z origin
  double thickness= slab.get_thickness();
  top_= 0.5*thickness;
  bottom_= -0.5*thickness;
  midZ_= 0;
  R_= slab.get_pore_radius(); // major radius
  rh_= slab.get_horizontal_minor_radius();
  rv_= slab.get_vertical_minor_radius();
  is_pore_radius_optimized_= slab.get_pore_radius_is_optimized();
}


//
inline double
SlabWithToroidalPorePairScore::evaluate_index
(Model *m,
 const ParticleIndexPair& pip,
 DerivativeAccumulator *da) const
{
  IMP_OBJECT_LOG;
  SlabWithToroidalPore slab(m, pip[0]);
  update_cached_slab_params(slab);
  IMP::core::XYZR d(m, pip[1]);
  if (!d.get_coordinates_are_optimized())
    return false;
  algebra::Sphere3D d_sphere( d.get_sphere() );
  algebra::Vector3D displacement;
  double score=get_sphere_penetration_depth(d.get_sphere(),
                               da ? &displacement : nullptr);
  IMP_LOG_TERSE("sphere " << d << " score " << score);
  if(da && score>0.0){
    algebra::Vector3D derivative_vector = -k_*displacement;
    IMP_LOG_TERSE(" derivative vector " << derivative_vector);
    d.add_to_derivatives(derivative_vector, *da);
    if(is_pore_radius_optimized_){
      double radial_displacement= displacement[0]*d_sphere[0]+displacement[1]*d_sphere[1]; // TODO: assumes slab origin at 0,0,0 - perhaps in the future extend to general case
      slab.add_to_pore_radius_derivative(k_*radial_displacement, *da);
    }

  }
  return score;
}

//
inline double
SlabWithToroidalPorePairScore::evaluate_indexes
    (Model *m,
     const ParticleIndexPairs &pips,
     DerivativeAccumulator *da,
     unsigned int lower_bound,
     unsigned int upper_bound) const
{
  IMP_LOG_TERSE("SlabWithToroidalPore singleton - evaluate indexes"
          << std::endl);
  if(upper_bound<lower_bound){
    return 0.0;
  }
  double ret(0.0);
    double radial_displacements(0.0); // sum of pore radius displacemnets
  algebra::Sphere3D const* spheres_table=
    m->access_spheres_data();
  algebra::Sphere3D* sphere_derivatives_table=
    m->access_sphere_derivatives_data();
  IMP::internal::BoolAttributeTableTraits::Container const&
    is_optimizable_table= m->access_optimizeds_data
    (core::XYZ::get_coordinate_key(0)); // x is indicator
  ParticleIndex slab_pi(pips[lower_bound][0]);
  SlabWithToroidalPore slab(m, slab_pi); // TODO: do this only in first round
  update_cached_slab_params(slab);

  // Evaluate and sum score and derivative for all particles:
  for (unsigned int i = lower_bound; i < upper_bound; ++i) {
    ParticleIndex pi( pips[i][1]);
    int pi_index=pi.get_index();
    // Check attributes have valid values:
    IMP_CHECK_CODE( {
        IMP_INTERNAL_CHECK(pips[i][0]==slab_pi,
                           "All particles are assumed to be evaluated against"
                           " the same slab");
        IMP::core::XYZR d(m, pi);
        algebra::Sphere3D s=spheres_table[pi_index];
        IMP_INTERNAL_CHECK(d.get_coordinates_are_optimized() == is_optimizable_table[pi],
                           "optimable table inconsistent with d.get_coordinates_are_optimized for particle " << d);
        IMP_INTERNAL_CHECK((d.get_coordinates() - s.get_center()).get_magnitude()<.001,
                           "Different coords for particle " << d << " *** "
                           << d.get_coordinates() << " vs. " << s.get_center());
        IMP_INTERNAL_CHECK(d.get_radius() == s.get_radius(),
                           "Different radii for particle " << d << " *** "
                           << d.get_radius() << " vs. " << s.get_radius());
      } ); // IMP_CHECK_CODE
    if(!is_optimizable_table[pi]) {
      continue;
    }
    algebra::Vector3D displacement;
    double cur_score = get_sphere_penetration_depth(spheres_table[pi_index],
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
      algebra::Sphere3D const& s=spheres_table[pi_index];
      radial_displacements+= displacement[0]*s[0] + displacement[1]*s[1]; // TODO: assumes slab origin at 0,0,0 - perhaps in דthe future extend to general case

    }
  }
  if(da && is_pore_radius_optimized_){
    slab.add_to_pore_radius_derivative(k_ * radial_displacements, *da);
  }
  return ret;
}

// sphere - sphere for which the distance is computed
// origin - ellipsoid origin (a point on the major circle of the torus)
double
SlabWithToroidalPorePairScore::
get_sphere_ellipsoid_penetration_depth
(algebra::Sphere3D const& sphere,
 algebra::Vector3D const& origin,
 algebra::Vector3D* out_translation) const
{
  const double eps= 1e-9;
  algebra::Vector3D v(sphere.get_center()-origin);
  double dXY2= v[0]*v[0]+v[1]*v[1];
  double dZ2= v[2]*v[2];
  double dv2= dXY2 + dZ2 + eps;
  // theta = atan(dXY/dZ) = asin(dXY/dv) = acos(dZ/dv)
  double sinTheta2= dXY2/dv2;
  double cosTheta2= dZ2/dv2;
  double cur_r= std::sqrt(rv_*rv_*cosTheta2 + rh_*rh_*sinTheta2);
  double dv=std::sqrt(dv2);
  double surface_distance= dv - sphere.get_radius() - cur_r;
  if(surface_distance>=0 ){
    if(out_translation){
      (*out_translation)=IMP::algebra::Vector3D(0,0,0);
    }
    return 0;
  }
  // Overlap:
  if(out_translation){
    (*out_translation)= v.get_unit_vector();
  }
  return -surface_distance; // penetration is negative the 'distance'


}

// return - distance to nearest surface point
// out_translation - a normalized vector pointing to the nearest surface point
double
SlabWithToroidalPorePairScore::get_sphere_penetration_depth
(algebra::Sphere3D sphere,
 algebra::Vector3D* out_translation) const
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
    if(out_translation){
      (*out_translation)=IMP::algebra::Vector3D(0,0,0);
    }
    return 0;
  }
  double d_xy2 = x*x + y*y; // distance from (0,0,z)
  bool is_closer_to_top= (sphere_dz_bottom > -sphere_dz_top);
  double abs_sphere_dz =
    is_closer_to_top ? -sphere_dz_top : sphere_dz_bottom; // d_z is the mianimal vertical change in z that brings s either fully above or fully below the slab
  if (d_xy2 > R_*R_) {
      // early about II - radially outside torus major circle of radius R
    if(out_translation){
      (*out_translation)=IMP::algebra::Vector3D(0,0,is_closer_to_top?+1:-1);
    }
    return abs_sphere_dz;
  }
  double d_xy = std::sqrt(d_xy2);
  // IV. In slab vertically + radially within torus major circle of radius R
  const double eps = 1e-9;
  double ret;
  if (d_xy>eps) {
    double scale= R_/d_xy;
    IMP::algebra::Vector3D origin(x*scale, y*scale, midZ_);
    ret= get_sphere_ellipsoid_penetration_depth
      (sphere, origin, out_translation );
  } else {
    IMP::algebra::Vector3D origin(R_, 0.0, midZ_);
    ret= get_sphere_ellipsoid_penetration_depth
      (sphere, origin, out_translation );
  }
  return ret;
}



IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_SLAB_WITH_TOROIDAL_PORE_PAIR_SCORE_H */
