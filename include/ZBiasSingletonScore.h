/**
 *  \file ZBiasSingletonScore.h
 *  \brief score that biases particles to go down the Z axis
 *
 *  Copyright 2007-8 Sali Lab. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_Z_BIAS_SINGLETON_SCORE_H
#define IMPNPCTRANSPORT_Z_BIAS_SINGLETON_SCORE_H

#include "npctransport_config.h"
#include "enums.h"
#include <IMP/SingletonScore.h>
#include <IMP/singleton_macros.h>
#include <IMP/algebra/Vector3D.h>
#include <limits>

IMPNPCTRANSPORT_BEGIN_NAMESPACE


//! Score that biases particles to go down the Z axis
class IMPNPCTRANSPORTEXPORT ZBiasSingletonScore
    : public SingletonScore {
 private:
  algebra::Vector3D v_deriv_; // deriv of force, which is of form (0,0,k)
  double max_r2_; // maximal square r (distance from origin on x,y plane)
 public:
  /**
     Exclude particles from the range of z coordinates [bottom_..top_]
     with repulsive force constant k

     @param k reuplsive force constant
     @param max_r maximal distance from origin on x,y plane (radius
                  relative to pore axis) in which force is applied.
                  if 0.0, no limiation
   */
 ZBiasSingletonScore
   ( double k,
     double max_r = HALF_SQRT_MAX_DOUBLE) // half for numerical margin error
   : max_r2_( max_r * max_r )
    {
      set_k(k);
    }


  /**
      returns the force constant for pulling to high z
      (negative = pull to low z)
  */
  double get_k() const { return v_deriv_[2]; }

  /**
      sets the force constant for pulling to high z
      (negative = pull to low z)

      @param k force contant
  */
  void set_k(double k) {
    v_deriv_ =  algebra::Vector3D(0, 0, k);
  }

  virtual double evaluate_index(Model *m, ParticleIndex pi,
                                DerivativeAccumulator *da) const IMP_OVERRIDE
  {
    core::XYZR d(m,pi);
    double r2 =  std::pow(d.get_x(), 2) + std::pow(d.get_y(), 2);
    if (r2 > max_r2_) {
      // nothing to add to derivatives( = 0);
      return 0.0; // TODO: this is not smooth around pore edge - problem
    }
    double const& k = v_deriv_[2];
    double score = k * d.get_z();
    if (da) {
      IMP_LOG(VERBOSE, "result in " << score
              << " and " << v_deriv_ << std::endl);
      d.add_to_derivatives(v_deriv_, *da);
    }
    return score;
  }

  virtual ModelObjectsTemp do_get_inputs(Model *m,
                                         const ParticleIndexes &pis) const
      IMP_OVERRIDE
  { return IMP::get_particles(m, pis);  }

  IMP_SINGLETON_SCORE_METHODS(ZBiasSingletonScore);
  IMP_OBJECT_METHODS(ZBiasSingletonScore);
};

IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_Z_BIAS_SINGLETON_SCORE_H */
