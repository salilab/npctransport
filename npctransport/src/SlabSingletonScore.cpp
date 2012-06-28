/**
 *  \file SeparateSingletonModifier.cpp
 *  \brief XXXX.
 *
 *  Copyright 2007-8 Sali Lab. All rights reserved.
 *
 */

#include "IMP/npctransport/SlabSingletonScore.h"
#include "IMP/core/XYZR.h"

IMPNPCTRANSPORT_BEGIN_NAMESPACE

// computes the distance and displacement vector of v
// from the surface of a z-axis aligned cylinder
//
// @return <distance, a vector pointing out>,
//         negative distance should mean v is inside the cylinder
std::pair<double, algebra::Vector3D>
SlabSingletonScore::get_displacement_vector(const algebra::Vector3D &v) const
{
  double R2= square(v[0]) + square(v[1]); // r^2 for [x,y] projection
  double H= v[2] - midZ_;  // height on z-axis from cyl origin
  //std::cout << h << " " << r2 << " for " << v << std::endl;
  if (R2 > square(radius_) || (v[2] <= top_ && v[2] >= bottom_))
    { // = either inside cylinder, or [x,y] outside cyl_radius
      double aH= std::abs(H); // absolute distance from cyl origin
      double R= std::sqrt(R2);
      double dR = R - radius_; // displacement on [x,y] direction
      if (dR + aH < .5 * height_) {
        //std::cout << "ring or tunnel" << std::endl;
        if (R2 < .00001) { // at origin
          return std::make_pair(radius_, algebra::Vector3D(0,0,1));
        } else {
          algebra::Vector3D rv(-v[0], -v[1], 0);
          return std::make_pair(radius_ - R, rv.get_unit_vector());
        }
      } else {
        //std::cout << "in or out of slab" << std::endl;
        if (H > 0) {
          return std::make_pair(v[2] - top_,    algebra::Vector3D(0,0,1));
        } else {
          return std::make_pair(bottom_ - v[2], algebra::Vector3D(0,0,-1));
        }
      }
    } else { // = outside cylinder && [x,y] within cyl_radius
      //std::cout << "channel" << std::endl;
      if (R2 < .00001) { // at origin
        //std::cout << "in center " << std::endl;
        if (H>0) {
          return std::make_pair(v[2] - top_, algebra::Vector3D(0,0,1));
        } else {
          return std::make_pair(bottom_ - v[2], algebra::Vector3D(0,0,-1));
        }
      }
      algebra::Vector3D rim
        =algebra::Vector3D(v[0], v[1], 0).get_unit_vector()*radius_;
      if (H > 0) {
        rim[2]= top_;
      } else {
        rim[2]= bottom_;
      }
      //std::cout << "rim is " << rim << std::endl;
      algebra::Vector3D diff= v-rim;
      return std::make_pair(diff.get_magnitude(), diff.get_unit_vector());
    }
}

SlabSingletonScore::SlabSingletonScore(double height, double radius, double k):
  height_(height), radius_(radius), k_(k),
   top_( height / 2.0 ), bottom_( -height / 2.0 ), midZ_( 0.0 )
{
}


double SlabSingletonScore::evaluate(Particle *p,
                                    DerivativeAccumulator *da) const {
  using algebra::Vector3D;
  IMP_OBJECT_LOG;
  core::XYZR d(p);
  if (!d.get_coordinates_are_optimized()) return false;
  // early abort if above or below slab
  if ( (d.get_z()-d.get_radius() > height_) ||
       (d.get_z()+d.get_radius() < -height_)) {
    return 0;
  }
  // early abort if [x,y] within cylinder radius
  if ((algebra::get_squared(d.get_x())+ algebra::get_squared(d.get_y()))
      < algebra::get_squared(radius_-d.get_radius())) {
    return 0;
  }
  std::pair<double, Vector3D> dp
    = get_displacement_vector(d.get_coordinates());
  IMP_LOG(VERBOSE, "At point " << d.get_coordinates()
          << " have distance " << dp.first << " and direction " << dp.second
          << std::endl);
  if (dp.first > d.get_radius())
    return 0;
  double distance = dp.first;
  Vector3D displacement = dp.second;
  IMP_INTERNAL_CHECK(std::abs(displacement.get_magnitude() -1) < .1,
                     "Not a unit vector");
  double score= k_*(d.get_radius() - distance);
  if (da) {
    Vector3D uv= -displacement;
    Vector3D dc=uv*k_;
    IMP_LOG(VERBOSE, "result in " << score << " and " << dc
            << std::endl);
    d.add_to_derivatives(dc, *da);
  }
  return score;
}
ParticlesTemp
SlabSingletonScore::get_input_particles(Particle*p) const {
  return ParticlesTemp(1,p);
}

ContainersTemp
SlabSingletonScore::get_input_containers(Particle*) const {
  return ContainersTemp();
}


void SlabSingletonScore::do_show(std::ostream &out) const {
  out << "SeparateSingletonModifier" << std::endl;
}

IMPNPCTRANSPORT_END_NAMESPACE
