/**
 *  \file SeparateSingletonScore.h
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
    for z=[-height..height]
    Returns 0 for all particles fully beyond [height] or below [-height],
    or fully within slab radius from the origin in the [X,Y] plabe
 */
class IMPNPCTRANSPORTEXPORT SlabSingletonScore: public SingletonScore
{
  double thickness_; // thichness of slab

  double radius_; // radius of slab cylinder

  double k_; // coefficient for violation of slab constraint

  double top_; // top of slab on z-axis

  double bottom_; // bottom of slab on x-axis

  double midZ_; // (top + bottom) / 2, for caching some calculations

 public:
  //! Get the individual particles from the passed SingletonContainer
  SlabSingletonScore(double thickness, double radius, double k);

  IMP_SINGLETON_SCORE(SlabSingletonScore);
  algebra::Vector3D
      get_displacement_direction(const algebra::Vector3D &v) const {
    return get_displacement_vector(v).second;
  }
  double get_displacement_magnitude(const algebra::Vector3D &v) const {
    return get_displacement_vector(v).first;
  }

 private:
  // computes the displacement of v from a z-axis aligned cylinder
  //
  // @return <distance, a vector pointing out>,
  //         negative distance means v is inside cylinder
  std::pair<double, algebra::Vector3D>
    get_displacement_vector(const algebra::Vector3D &v) const;
};


IMPNPCTRANSPORT_END_NAMESPACE

#endif  /* IMPNPCTRANSPORT_SLAB_SINGLETON_SCORE_H */
