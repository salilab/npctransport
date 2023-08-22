/**
 *  \file SlabWithCylindricalPorePairScore.h
 *  \brief XXXXXXXXXXXXXX
 *
 *  Copyright 2007-2022 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_SLAB_PAIR_SCORE_H
#define IMPNPCTRANSPORT_SLAB_PAIR_SCORE_H

#include "npctransport_config.h"
#include "SlabWithCylindricalPore.h"
#include <IMP/check_macros.h>
#include <IMP/PairScore.h>
#include <IMP/pair_macros.h>
#include "IMP/core/XYZR.h"

IMPNPCTRANSPORT_BEGIN_NAMESPACE

//! XXXX
/** An origin centered slab with a pore in the vertical direction,
    for z = [-0.5*thickness_...0.5*thickness_]

    The score evaluates to 0 for all particles fully beyond z range
    or fully within slab radius from the origin in the [X,Y] plane
    For particles that penetrate the slab, the score increases linearly
    with the magnitude of penetration with a slope k in units of kcal/mol/A.
    Conseuqently, the gradient is constant and it is oriented so as to repulse
    the particle towards the nearest point on the slab surface.
 */
class IMPNPCTRANSPORTEXPORT
SlabWithCylindricalPorePairScore : public PairScore {
 private:
  double k_;  // coefficient for violation of a slab restraint in kcal/mol/A

  // cache variables (therefore, all are mutable, as they are only used for performance purposes)
  mutable double thickness_;  // thickness of slab
  mutable double pore_radius_;  // radius of slab cylinder
  mutable double top_;  // top of slab on z-axis
  mutable double bottom_;  // bottom of slab on x-axis
  mutable double midZ_;  // (top + bottom) / 2, for caching some calculations
  mutable bool is_pore_radius_optimized_;

 public:
  //! Constructs a slab pair score that acts on a
  //! SlabWithCylindricalPore particle and a diffusing particle.
  //! The slab applies a repulsive force constant k specified in units
  //! of kcal/mol/A (linear potential)
  SlabWithCylindricalPorePairScore(double k);

  //! returns the direction vector for the displacement of point v relative to the slab surface
    algebra::Vector3D get_displacement_direction
    (SlabWithCylindricalPore const& slab, const algebra::Vector3D &v) const;

  //! returns the displacement magnitude of point v relative to the slab surface
  double get_displacement_magnitude
    (SlabWithCylindricalPore const&slab, const algebra::Vector3D &v) const;

 public:
  //! evaluate score for particle pair pip in model m
  /** evaluate score for particle pair pip in model m, where the first particle
      is assumed to be a cylindrical slab. If da is not null,
      use it to accumulate derivatives in the model.
  */
  virtual double evaluate_index
    (Model *m,
     const ParticleIndexPair& pip,
     DerivativeAccumulator *da) const override;

  virtual ModelObjectsTemp do_get_inputs(Model *m,
                                         const ParticleIndexes &pis) const
      override;

  //! evaluate score for particles pis[lower_bound..upper_bound] in
  //! model m. If da is not null, use it to accumulate derivatives
  //! in model.
  virtual double evaluate_indexes
    (Model *m,
     const ParticleIndexPairs &pips,
     DerivativeAccumulator *da,
     unsigned int lower_bound,
     unsigned int upper_bound) const override final;

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
      const ParticleIndexPairs &pips,
      DerivativeAccumulator *da,
      double max,
      unsigned int lower_bound,
      unsigned int upper_bound) const override
  {
    double ret = 0;
    for (unsigned int i = lower_bound; i < upper_bound; ++i) {
      ret += evaluate_if_good_index(m, pips[i], da, max - ret);
      if (ret > max) return std::numeric_limits<double>::max();
    }
    return ret;
  }

  //  IMP_PAIR_SCORE_METHODS(SlabWithCylindricalPorePairScore);
  IMP_OBJECT_METHODS(SlabWithCylindricalPorePairScore);

 private:
  // evaluate slab for specified sphere, based on most recent cached
  // slab parameters (/see update_cached_slab_params()).
  // Return 0 if ball does not penetrate slab.
  //
  // @param s the sphere to evaluate
  // @param out_displacement if not null and the returned score is positive, *out_displacement
  //                         is used to store the computed displacement vector from
  //                         the surface of the z-axis aligned cylinder to
  //                         the center of s. Ignore if score is zero.
  //
  inline double evaluate_sphere
    (algebra::Sphere3D s,
     algebra::Vector3D* out_displacement) const;


  // computes the displacement from the surface of the z-axis
  // aligned cylinder to v, based on the most recent cached
  // slab parameters (/see update_cached_slab_params()).
  //
  // @return <distance, a vector pointing out>,
  //         negative distance means v is inside cylinder
  std::pair<double, algebra::Vector3D> get_displacement_vector(
      const algebra::Vector3D &v) const;

  // update internal variables holding slab params for fast access
  // based on decorated particle slab
  void update_cached_slab_params
    (SlabWithCylindricalPore slab) const;
};

inline void
SlabWithCylindricalPorePairScore::update_cached_slab_params
(SlabWithCylindricalPore slab) const
{
  //TODO: support slabs with non-zero x,y,z origin
  thickness_= slab.get_thickness();
  top_= 0.5*thickness_;
  bottom_= -0.5*thickness_;
  midZ_= 0;
  pore_radius_= slab.get_pore_radius();
  is_pore_radius_optimized_= slab.get_pore_radius_is_optimized();
}

//
inline double
SlabWithCylindricalPorePairScore::evaluate_index
(Model *m,
 const ParticleIndexPair& pip,
 DerivativeAccumulator *da) const
{
  IMP_OBJECT_LOG;
  IMP_USAGE_CHECK(SlabWithCylindricalPore::get_is_setup(m, pip[0]),
                  "pip[0] is not a SlabWithCylindricalPore in evaluate_index()");
  SlabWithCylindricalPore slab(m, pip[0]);
  update_cached_slab_params(slab);
  IMP::core::XYZR d(m, pip[1]);
  algebra::Sphere3D d_sphere( d.get_sphere() );
  if (!d.get_coordinates_are_optimized())
    return false;
  algebra::Vector3D displacement; //  a unit displacement vector - output of evaluate sphere
  double score=evaluate_sphere(d_sphere,
                               da ? &displacement : nullptr);
  if(da && score>0.0){
    algebra::Vector3D derivative_vector = -k_*displacement;
    IMP_LOG(PROGRESS, "result in " << score << " and " << derivative_vector << std::endl);
    d.add_to_derivatives(derivative_vector, *da);
    if(is_pore_radius_optimized_){
      // TODO: assume that the direction of a positive radial displacement vector is opposite to the sphere x,y vector - is this always true?
      double radial_displacement_magnitude= // magnitude of the displacement vector projected on the x,y plane
        std::sqrt(displacement[0]*displacement[0]+displacement[1]*displacement[1]); // TODO: currently we assume slab origin at 0,0,0 - perhaps in the future extend to general case
      slab.add_to_pore_radius_derivative(-k_*radial_displacement_magnitude, *da);
    }
  }
  return score;
}

//
inline double
SlabWithCylindricalPorePairScore::evaluate_indexes
(Model *m,
 const ParticleIndexPairs &pips,
 DerivativeAccumulator *da,
 unsigned int lower_bound,
 unsigned int upper_bound) const
{
  if(upper_bound<lower_bound){
    return 0.0;
  }
  double ret(0.0);
  double radial_displacements_magnitude(0.0); // sum of pore radius displacemnets
  algebra::Sphere3D const* spheres_table=
    m->access_spheres_data();
  algebra::Sphere3D* sphere_derivatives_table=
    m->access_sphere_derivatives_data();
  IMP::internal::BoolAttributeTableTraits::Container const& is_optimizable_table=
    m->access_optimizeds_data(core::XYZ::get_coordinate_key(0)); // use only x coordinate as indicator
  ParticleIndex slab_pi(pips[lower_bound][0]);
  SlabWithCylindricalPore slab(m, slab_pi); // TODO: do this only in first round
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
    double cur_score = evaluate_sphere(spheres_table[pi_index],
                                        da ? &displacement : nullptr);
    ret+= cur_score;
    if(cur_score>0.0 && da) {
      algebra::Vector3D derivative_vector = -k_*displacement;
      // accumulate derivatives directly for speed
      for(unsigned int j=0; j<3; j++) {
        sphere_derivatives_table[pi_index][j] += (*da)(derivative_vector[j]);
      }
      // TODO: assume that the direction of a positive radial displacement vector is opposite to the sphere x,y vector - is this always true?
      radial_displacements_magnitude+=
        std::sqrt(displacement[0]*displacement[0] + displacement[1]*displacement[1]); // TODO: assumes slab origin at 0,0,0 - perhaps in דthe future extend to general case
    }
  }
  if(da && is_pore_radius_optimized_){
    slab.add_to_pore_radius_derivative(-k_ * radial_displacements_magnitude, *da);
  }
  return ret;
}

//
double
SlabWithCylindricalPorePairScore::evaluate_sphere
(algebra::Sphere3D s,
 algebra::Vector3D* out_displacement) const
{
  IMP_OBJECT_LOG;
  IMP_LOG(VERBOSE, "evaluate_sphere " << s << std::endl);
  double const x= s[0];
  double const y= s[1];
  double const z= s[2];
  double const sr= s.get_radius();
  // early abort if above or below slab
  if ((z-sr > top_) || (z+sr < bottom_)) {
    return 0;
  }
  double const x2= x*x;
  double const y2= y*y;
  double const R= pore_radius_-sr;
  double const R2= R*R;
  // early abort if [x,y] within cylinder perimeter
  if (x2+y2 < R2) {
    return 0;
  }
  std::pair<double, algebra::Vector3D> dp = get_displacement_vector(s.get_center());
  IMP_LOG(PROGRESS,
          "At point " << s.get_center() << " have distance " << dp.first
          << " and direction " << dp.second << std::endl);
  double const distance = dp.first;
  if (distance > sr) {
    return 0;
  }
  double const score= k_ * (sr - distance); // must be positive if distance <= sr
  if(out_displacement){
    *out_displacement= dp.second;
    IMP_INTERNAL_CHECK(std::abs(out_displacement->get_magnitude() - 1) < .1,
                       "Not a unit vector");
  }
  return score;
}

// computes the distance and a unit displacement vector of v
// from the surface of a z-axis aligned cylinder
//
// @return <distance, a unit vector pointing outwards>,
//         negative distance should mean v is inside the cylinder
inline std::pair<double, algebra::Vector3D>
SlabWithCylindricalPorePairScore::get_displacement_vector(const algebra::Vector3D &v) const {
  double dXY2 = square(v[0]) + square(v[1]);  // r^2 for [x,y] projection
  double dZ = v[2] - midZ_;  // thickness on z-axis from cyl origin
  IMP_LOG_PROGRESS( dZ << " " << dXY2 << " for " << v << std::endl);
  if (dXY2 > square(pore_radius_) ||
      (v[2] <= top_ && v[2] >= bottom_)) {  //  inside vertical slab boundaries and pore perimeter on x,y plane OR outside pore perimeter in any vertical position
    double abs_dZ = std::abs(dZ);
    double abs_dXY = std::sqrt(dXY2);
    double dR = abs_dXY - pore_radius_;  // displacement on [x,y] direction (positive = outside pore perimeter on x,y plane)
    if (dR + abs_dZ < .5 * thickness_) { // in a cones from (0,0,.5thickness) to (0,.5thicness,0)
      IMP_LOG_PROGRESS("ring or pore" << std::endl);
      if (dXY2 < .00001) {  // at origin
        return std::make_pair(pore_radius_, algebra::Vector3D(0, 0, 1));
      } else {
        algebra::Vector3D rv(-v[0], -v[1], 0);
        return std::make_pair(pore_radius_ - abs_dXY, rv.get_unit_vector());
      }
    } else { // possibly in the pore but out of the 'double-cone diamond'
      IMP_LOG_PROGRESS("in or out of slab" << std::endl);
      if (dZ > 0) {
        return std::make_pair(v[2] - top_, algebra::Vector3D(0, 0, 1));
      } else {
        return std::make_pair(bottom_ - v[2], algebra::Vector3D(0, 0, -1));
      }
    }
  } else {  // = outside slab boundaries AND insider pore perimeter on x,y plane
    IMP_LOG_PROGRESS("channel" << std::endl);
    if (dXY2 < .00001) {  // at central axis
      IMP_LOG_PROGRESS("in center " << std::endl);
      if (dZ > 0) {
        return std::make_pair(v[2] - top_, algebra::Vector3D(0, 0, 1));
      } else {
        return std::make_pair(bottom_ - v[2], algebra::Vector3D(0, 0, -1));
      }
    }
    algebra::Vector3D rim =
        algebra::Vector3D(v[0], v[1], 0).get_unit_vector() * pore_radius_;
    if (dZ > 0) {
      rim[2] = top_;
    } else {
      rim[2] = bottom_;
    }
    IMP_LOG_PROGRESS( "rim is " << rim << std::endl);
    algebra::Vector3D diff = v - rim;
    return std::make_pair(diff.get_magnitude(), diff.get_unit_vector());
  }
}




IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_SLAB_PAIR_SCORE_H */
