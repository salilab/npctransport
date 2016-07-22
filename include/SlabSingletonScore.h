/**
 *  \file SlabSingletonScore.h
 *  \brief XXXXXXXXXXXXXX
 *
 *  Copyright 2007-8 Sali Lab. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_SLAB_SINGLETON_SCORE_H
#define IMPNPCTRANSPORT_SLAB_SINGLETON_SCORE_H

#include "npctransport_config.h"
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

  virtual double evaluate_index(Model *m, ParticleIndex p,
                                DerivativeAccumulator *da) const IMP_OVERRIDE;
  virtual ModelObjectsTemp do_get_inputs(Model *m,
                                         const ParticleIndexes &pis) const
      IMP_OVERRIDE;

  /* virtual double evaluate_indexes(Model *m, const ParticleIndexPairs &p, */
  /*                                 DerivativeAccumulator *da, */
  /*                                 unsigned int lower_bound, */
  /*                                 unsigned int upper_bound) const IMP_FINAL; */

  /* double evaluate_if_good_indexes */
  /*   ( Model *m, const ParticleIndexes &p, DerivativeAccumulator *da,            */
  /*     double max, unsigned int lower_bound, unsigned int upper_bound) const */
  /* {   */
  /*   double ret = 0;                                                             */
  /*   for (unsigned int i = lower_bound; i < upper_bound; ++i) {                  */
  /*     ret += evaluate_if_good_index(m, p[i], da, max - ret);                    */
  /*     if (ret > max) return std::numeric_limits<double>::max();                 */
  /*   }                                                                           */
  /*   return ret;                                                                 */
  /* } */

  IMP_SINGLETON_SCORE_METHODS(SlabSingletonScore);
  IMP_OBJECT_METHODS(SlabSingletonScore);

 private:
  // computes the displacement of v from a z-axis aligned cylinder
  //
  // @return <distance, a vector pointing out>,
  //         negative distance means v is inside cylinder
  std::pair<double, algebra::Vector3D> get_displacement_vector(
      const algebra::Vector3D &v) const;
};

inline double
SlabSingletonScore::evaluate_index(Model *m, const ParticleIndex pi,
                                   DerivativeAccumulator *da) const {
  using algebra::Vector3D;
  IMP_OBJECT_LOG;
  IMP::core::XYZR d(m, pi);
  if (!d.get_coordinates_are_optimized()) return false;
  algebra::Sphere3D const s=d.get_sphere();
  double z=s[2];
  double sr=s.get_radius();
  // early abort if above or below slab
  if ((z-sr > top_) || (z+sr < bottom_)) {
    return 0;
  }
  double x2=algebra::get_squared(s[0]);
  double y2=algebra::get_squared(s[1]);
  double r2=algebra::get_squared(radius_-sr);
  // early abort if [x,y] within cylinder radius
  if (x2+y2 < r2) {
    return 0;
  }
  std::pair<double, Vector3D> dp = get_displacement_vector(s.get_center());
  IMP_LOG(VERBOSE,
          "At point " << d.get_coordinates() << " have distance " << dp.first
          << " and direction " << dp.second << std::endl);
  double distance = dp.first;
  Vector3D displacement = dp.second;
  if (distance > sr) return 0;
  IMP_INTERNAL_CHECK(std::abs(displacement.get_magnitude() - 1) < .1,
                     "Not a unit vector");
  double score = k_ * (sr - distance);
  if (da) {
    Vector3D uv = -displacement;
    Vector3D dc = uv * k_;
    IMP_LOG(VERBOSE, "result in " << score << " and " << dc << std::endl);
    d.add_to_derivatives(dc, *da);
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
  // std::cout << h << " " << r2 << " for " << v << std::endl;
  if (R2 > square(radius_) ||
      (v[2] <= top_ && v[2] >= bottom_)) {  // = either inside cylinder, or
                                            // [x,y] outside cyl_radius
    double aH = std::abs(H);  // absolute distance from cyl origin
    double R = std::sqrt(R2);
    double dR = R - radius_;  // displacement on [x,y] direction
    if (dR + aH < .5 * thickness_) {
      // std::cout << "ring or tunnel" << std::endl;
      if (R2 < .00001) {  // at origin
        return std::make_pair(radius_, algebra::Vector3D(0, 0, 1));
      } else {
        algebra::Vector3D rv(-v[0], -v[1], 0);
        return std::make_pair(radius_ - R, rv.get_unit_vector());
      }
    } else {
      // std::cout << "in or out of slab" << std::endl;
      if (H > 0) {
        return std::make_pair(v[2] - top_, algebra::Vector3D(0, 0, 1));
      } else {
        return std::make_pair(bottom_ - v[2], algebra::Vector3D(0, 0, -1));
      }
    }
  } else {  // = outside cylinder && [x,y] within cyl_radius
    // std::cout << "channel" << std::endl;
    if (R2 < .00001) {  // at origin
      // std::cout << "in center " << std::endl;
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
    // std::cout << "rim is " << rim << std::endl;
    algebra::Vector3D diff = v - rim;
    return std::make_pair(diff.get_magnitude(), diff.get_unit_vector());
  }
}




IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_SLAB_SINGLETON_SCORE_H */
