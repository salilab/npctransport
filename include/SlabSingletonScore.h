/**
 *  \file SlabSingletonScore.h
 *  \brief XXXXXXXXXXXXXX
 *
 *  Copyright 2007-8 Sali Lab. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_SLAB_SINGLETON_SCORE_H
#define IMPNPCTRANSPORT_SLAB_SINGLETON_SCORE_H

#include "npctransport_config.h"
#include <IMP/check_macros.h>
#include <IMP/SingletonScore.h>
#include <IMP/singleton_macros.h>
#include "IMP/core/XYZR.h"

IMPNPCTRANSPORT_BEGIN_NAMESPACE

//! XXXX
/** An origin centered slab with a tunnel in the vertical direction,
    for z = [-0.5*thickness_...0.5*thickness_]
    Returns 0 for all particles fully beyond z range
    or fully within slab radius from the origin in the [X,Y] plane
    // TODO: verify documentation
 */
class IMPNPCTRANSPORTEXPORT SlabSingletonScore : public SingletonScore {
  double thickness_;  // thichness of slab

  double radius_;  // radius of slab cylinder

  double k_;  // coefficient for violation of slab constraint

  double top_;  // top of slab on z-axis

  double bottom_;  // bottom of slab on x-axis

  double midZ_;  // (top + bottom) / 2, for caching some calculations

 public:
  //! Get the individual particles from the passed SingletonContainer
  SlabSingletonScore(double thickness, double radius, double k);

  algebra::Vector3D get_displacement_direction(
      const algebra::Vector3D &v) const {
    return get_displacement_vector(v).second;
  }
  double get_displacement_magnitude(const algebra::Vector3D &v) const {
    return get_displacement_vector(v).first;
  }

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

  //  IMP_SINGLETON_SCORE_METHODS(SlabSingletonScore);
  IMP_OBJECT_METHODS(SlabSingletonScore);

 private:
  // evaluate slab for specified sphere. Return 0 if ball
  // does not penetrate slab
  //
  // @param s the sphere to evaluate
  // @param out_displacement if not null and the returned score is positive, *out_displacement
  //                         is used to store the computed displacement vector from
  //                         the surface of the z-axis aligned cylinder to
  //                         the center of s. Ignore if score is zero.
  double evaluate_sphere
    (algebra::Sphere3D s,
     algebra::Vector3D* out_displacement) const;


  // computes the displacement from the surface of the z-axis
  // aligned cylinder to v
  //
  // @return <distance, a vector pointing out>,
  //         negative distance means v is inside cylinder
  std::pair<double, algebra::Vector3D> get_displacement_vector(
      const algebra::Vector3D &v) const;
};

inline double
SlabSingletonScore::evaluate_index
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

inline double
SlabSingletonScore::evaluate_indexes(Model *m, const ParticleIndexes &pis,
                                DerivativeAccumulator *da,
                                unsigned int lower_bound,
                                unsigned int upper_bound) const
{
  double ret = 0;
  // Direct access to pertinent attributes:
  algebra::Sphere3D const* spheres_table=
    m->access_spheres_data();
  algebra::Sphere3D* sphere_derivatives_table=
    m->access_sphere_derivatives_data();
  IMP::internal::BoolAttributeTableTraits::Container const& is_optimizable_table=
    m->access_optimizeds_data(core::XYZ::get_coordinate_key(0)); // use only x coordinate as indicator
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


inline double
SlabSingletonScore::evaluate_sphere(algebra::Sphere3D s,
                                    algebra::Vector3D* out_displacement) const {
  double const x=s[0];
  double const y=s[1];
  double const z=s[2];
  double const sr=s.get_radius();
  // early abort if above or below slab
  if ((z-sr > top_) || (z+sr < bottom_)) {
    return 0;
  }
  double const x2=x*x;
  double const y2=y*y;
  double const R=radius_-sr;
  double const R2=R*R;
  // early abort if [x,y] within cylinder radius
  if (x2+y2 < R2) {
    return 0;
  }
  std::pair<double, algebra::Vector3D> dp = get_displacement_vector(s.get_center());
  IMP_LOG(VERBOSE,
          "At point " << s.get_center() << " have distance " << dp.first
          << " and direction " << dp.second << std::endl);
  double const distance = dp.first;
  if (distance > sr) {
    return 0;
  }
  double const score = k_ * (sr - distance); // must be positive score now
  if(out_displacement){
    *out_displacement = dp.second;
    IMP_INTERNAL_CHECK(std::abs(out_displacement->get_magnitude() - 1) < .1,
                       "Not a unit vector");
  }
  return score;
}



// computes the distance and displacement vector of v
// from the surface of a z-axis aligned cylinder
//
// @return <distance, a vector pointing out>,
//         negative distance should mean v is inside the cylinder
inline std::pair<double, algebra::Vector3D>
SlabSingletonScore::get_displacement_vector(const algebra::Vector3D &v) const {
  double R2 = square(v[0]) + square(v[1]);  // r^2 for [x,y] projection
  double H = v[2] - midZ_;  // thickness on z-axis from cyl origin
  IMP_LOG_PROGRESS( H << " " << R2 << " for " << v << std::endl);
  if (R2 > square(radius_) ||
      (v[2] <= top_ && v[2] >= bottom_)) {  // = either inside cylinder, or
                                            // [x,y] outside cyl_radius
    double aH = std::abs(H);  // absolute distance from cyl origin
    double R = std::sqrt(R2);
    double dR = R - radius_;  // displacement on [x,y] direction
    if (dR + aH < .5 * thickness_) {
      std::cout << "ring or tunnel" << std::endl;
      if (R2 < .00001) {  // at origin
        return std::make_pair(radius_, algebra::Vector3D(0, 0, 1));
      } else {
        algebra::Vector3D rv(-v[0], -v[1], 0);
        return std::make_pair(radius_ - R, rv.get_unit_vector());
      }
    } else {
      IMP_LOG_PROGRESS("in or out of slab" << std::endl);
      if (H > 0) {
        return std::make_pair(v[2] - top_, algebra::Vector3D(0, 0, 1));
      } else {
        return std::make_pair(bottom_ - v[2], algebra::Vector3D(0, 0, -1));
      }
    }
  } else {  // = outside cylinder && [x,y] within cyl_radius
    IMP_LOG_PROGRESS("channel" << std::endl);
    if (R2 < .00001) {  // at origin
      IMP_LOG_PROGRESS("in center " << std::endl);
      if (H > 0) {
        return std::make_pair(v[2] - top_, algebra::Vector3D(0, 0, 1));
      } else {
        return std::make_pair(bottom_ - v[2], algebra::Vector3D(0, 0, -1));
      }
    }
    algebra::Vector3D rim =
        algebra::Vector3D(v[0], v[1], 0).get_unit_vector() * radius_;
    if (H > 0) {
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

#endif /* IMPNPCTRANSPORT_SLAB_SINGLETON_SCORE_H */
