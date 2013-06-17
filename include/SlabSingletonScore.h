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

IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_SLAB_SINGLETON_SCORE_H */
