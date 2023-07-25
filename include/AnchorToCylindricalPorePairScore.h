/**
 *  \file AnchorToCylindricalPorePairScore.h
 *  \brief A chore for anchoring a bead to a pore
 *
 *  Copyright 2007-2022 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_ANCHOR_TO_CYLINDRICAL_PORE_PAIR_SCORE_H
#define IMPNPCTRANSPORT_ANCHOR_TO_CYLINDRICAL_PORE_PAIR_SCORE_H

#include "npctransport_config.h"
#include "SlabWithCylindricalPore.h"
#include <IMP/check_macros.h>
#include <IMP/PairScore.h>
#include <IMP/pair_macros.h>
#include <IMP/core/XYZR.h>
#include <IMP/score_functor/Harmonic.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

/**
   a score for anchoring a bead to  a reference point defined relatively to the surface of a cylindrical pore. The reference point is at a fixed angle about the pore central axis, radial distance from the pore edge, and z coordinate.
*/
class IMPNPCTRANSPORTEXPORT
AnchorToCylidnricalPorePairScore : public PairScore
{
 private:
  algebra::Vector2D normalized_xy_; // normalized direction of equilibrium point in x,y plane
  Float pore_radial_d_; // equilibrium distance of anchor point from pore surface in x,y plane
  score_functor::Harmonic ds_; // distance score
  mutable algebra::Vector3D reference_point_; // current value of equilibrium reference point

 public:
  /** Initialize a harmonic distance score between an anchor bead
      and an equilibrium reference point that is evaluated according to the specified
      equilibrium values, relative to the cylindrical pore at the time of score
      evaluation.

      @param rot_angle equilibrium counter-clockwise angle of the anchor bead
                       in the x,y plane relative to the x axis, in Radians
      @param radial_d equilibrium radial distance from the pore surface (positive = inside pore; negative = outside)
      @param z equilibrium z coordinate
      @param k the harmonic force coefficient between the bead and the equilibrium
               reference point
  */
  AnchorToCylidnricalPorePairScore(Float rot_angle,
                                   Float pore_radial_d,
                                   Float z,
                                   Float k);

  /** initialized the score using current anchor_bead orientation
      relative to scp for initializing the values of the equilibrium
      rotation angle of the anchor about the pore central axis,
      its radial distance from the pore surface, and its z coordinate

      @param scp The porus slab
      @param initial_anchor_point The coordinates that are used as reference
                 for initializing the anchor equilibrium values relative to the slab
      @param k   The harmonic force coefficient between the bead and the equilibrium
                 reference point
  */
  AnchorToCylidnricalPorePairScore(SlabWithCylindricalPore scp,
                                   algebra::Vector3D initial_anchor_point,
                                   Float k);

 public:
  //! Evaluate score for particle pair pip in model m
  /**
     Evaluate score for particle pair pip in model m, where pip[0] is
     a SlabWithCylindricalPore and pip[1] is the anchored bead.
     The score is a harmonic distance restraint between the anchored bead
     and the equilibrium reference point relative to the porous slab
     (see constructors for details).

     If da is not null, use it to accumulate derivatives in model.
  */
  virtual double evaluate_index
    (Model *m,
     const ParticleIndexPair& pip,
     DerivativeAccumulator *da) const override;

  virtual ModelObjectsTemp do_get_inputs(Model *m,
                                         const ParticleIndexes &pis) const
      override;

  IMP_PAIR_SCORE_METHODS(AnchorToCylidnricalPorePairScore);
  IMP_OBJECT_METHODS(AnchorToCylidnricalPorePairScore);

 private:

  // update internal variable holding current reference point for fast access
  // based on decorated particle slab
  // @param pr current pore radius
  void update_reference_point_for_pore_radius
    (Float pr) const;
};

inline void
AnchorToCylidnricalPorePairScore
::update_reference_point_for_pore_radius
( Float pr) const
{
  Float r= pr - pore_radial_d_;
  reference_point_[0]= normalized_xy_[0]*r;
  reference_point_[1]= normalized_xy_[1]*r;
  }

//
inline double
AnchorToCylidnricalPorePairScore
::evaluate_index
( Model *m,
  const ParticleIndexPair& pip,
  DerivativeAccumulator *da ) const
{
  IMP_USAGE_CHECK(SlabWithCylindricalPore::get_is_setup(m, pip[0]),
                  "pip[0] is not a SlabWithCylindricalPore in evaluate_index()");
  SlabWithCylindricalPore slab(m, pip[0]);
  update_reference_point_for_pore_radius(slab.get_pore_radius());
  IMP::core::XYZ xyz(m, pip[1]);
  if(!xyz.get_coordinates_are_optimized()){
    return false;
  }
  algebra::Vector3D rp_to_xyz =
    reference_point_ - xyz.get_coordinates();
  double sq = rp_to_xyz.get_squared_magnitude();
  double dist = std::sqrt(sq);
  if (da) {
    std::pair<double, double> sp =
      ds_.get_score_and_derivative(m, pip, dist);
    static const double MIN_DISTANCE = .00001;
    algebra::Vector3D rp_to_xyz_norm;
    if (dist > MIN_DISTANCE) {
      rp_to_xyz_norm = rp_to_xyz / dist;
    } else {
      rp_to_xyz_norm = algebra::get_zero_vector_d<3>();
    }
    algebra::Vector3D deriv_on_rp(rp_to_xyz_norm * sp.second);
    if(slab.get_pore_radius_is_optimized()){
      // transform derivative of slab reference position to derivative of pore radius
      double deriv_on_pore_radius_squared=
        deriv_on_rp[0]*deriv_on_rp[0] + deriv_on_rp[1]*deriv_on_rp[1];
      slab.add_to_pore_radius_derivative(sqrt(deriv_on_pore_radius_squared), *da);
    }
    m->add_to_coordinate_derivatives(pip[1], -deriv_on_rp, *da);
    return sp.first;
  } else {
    return ds_.get_score(m, pip, dist);
  }
}




IMPNPCTRANSPORT_END_NAMESPACE

#endif /*  IMPNPCTRANSPORT_ANCHOR_TO_CYLINDRICAL_PORE_PAIR_SCORE_H */
