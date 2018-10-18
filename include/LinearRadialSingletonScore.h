/**
 *  \file RadialSingletongScore.h
 *  \brief score that pulls particles towards the periphery
 *
 *  Copyright 2007-8 Sali Lab. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_RADIAL_SINGLETON_SCORE_H
#define IMPNPCTRANSPORT_RADIAL_SINGLETON_SCORE_H

#include "npctransport_config.h"
#include "enums.h"
#include <IMP/SingletonScore.h>
#include <IMP/singleton_macros.h>
#include <IMP/algebra/Vector3D.h>
#include <limits>

IMPNPCTRANSPORT_BEGIN_NAMESPACE


/**
   A singleton score that increases or decreases linearly
   as particles get farther away from the center of some
   sphere
*/
class IMPNPCTRANSPORTEXPORT
LinearRadialSingletonScore
: public IMP::SingletonScore
{
 private:
  //! The score increases or decreases radially from this center
  algebra::Vector3D center_;

  //! the force constant (in kcal/mol/A for a linear score)
  double k_;
 public:
  /**
     Push particles radially away or towards the center of some sphere

     @param center the sphere center of the radial score
     @param k the magnitude of the radial force vector in kcal/mol/A.
            It is positive for a divergent force and negative for a
            convergent (gravitation-like) force.  In other words, k_
            is the negative gradient of the score as the radius
            increases.
   */
  LinearRadialSingletonScore
   ( algebra::Vector3D center,
     double k)
    : center_( center),  k_( k )
    {
    }


  /**
      returns the force constant for pushing particles away
      from the center in kcal/mol/A (negative = pull)
  */
  double get_k() const {
    return k_;
  }

  /**
     sets the force constant for pushing particles away
     from the center in kcal/mol/A (negative=pull)

     @param k force contant
  */
  void set_k(double k) {
    k_= k;
  }

  virtual double evaluate_index(Model *m, ParticleIndex pi,
                                DerivativeAccumulator *da) const IMP_OVERRIDE
  {
    core::XYZ xyz(m,pi);
    algebra::Vector3D dxyz= xyz.get_coordinates() - center_;
    double dxyz_magnitude=
      algebra::get_magnitude_and_normalize_in_place( dxyz );
    // TODO: BR did not implement the dependence on the number of bound patches, so here is a good place
    //       to read it and scale by it
    double score = -k_ * dxyz_magnitude; // score decreases (improves) radially
    if (da) {
      algebra::Vector3D& dxyz_normalized= dxyz; // it is now a normalized version of itself
      algebra::Vector3D deriv= -k_ * dxyz_normalized;
      IMP_LOG(VERBOSE, "score " << score
              << " and derivative " << deriv << std::endl);
      xyz.add_to_derivatives(deriv, *da);
    }
    return score;
  }

  virtual ModelObjectsTemp do_get_inputs(Model *m,
                                         const ParticleIndexes &pis) const
      IMP_OVERRIDE
  { return IMP::get_particles(m, pis);  }

  IMP_SINGLETON_SCORE_METHODS(LinearRadialSingletonScore);
  IMP_OBJECT_METHODS(LinearRadialSingletonScore);
};

IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_LINEAR_RADIAL_SINGLETON_SCORE_H */
