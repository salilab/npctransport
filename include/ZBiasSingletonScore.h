/**
 *  \file ZBiasSingletonScore.h
 *  \brief score that biases particles to go down the Z axis
 *
 *  Copyright 2007-2022 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_Z_BIAS_SINGLETON_SCORE_H
#define IMPNPCTRANSPORT_Z_BIAS_SINGLETON_SCORE_H

#include "npctransport_config.h"
#include "enums.h"
#include <IMP/SingletonScore.h>
#include <IMP/singleton_macros.h>
#include <IMP/algebra/Vector3D.h>
#include <limits>
#include <utility>

IMPNPCTRANSPORT_BEGIN_NAMESPACE


//! Score that biases particles to go down the Z axis
class IMPNPCTRANSPORTEXPORT ZBiasSingletonScore
    : public SingletonScore {
 private:
  double k_; // force constant (positive for pulling towards z_, negative for pushing away from z_)
  double max_r2_; // maximal square r (distance from z axis)
  double z_; // z to bias towards/away from
 public:
  /**
     Exclude particles from the range of z coordinates [bottom_..top_]
     with repulsive force constant k

     @param k force constant (positive for pulling towards z, negative for pushing away from z)
     @param max_r maximal distance from z axis (radius
                  relative to pore axis) in which force is applied.
                  if 0.0, no limitation
     @param z z of xy plane to bias towards / against (relevant only if k is non-zero)
   */
 ZBiasSingletonScore
   ( double k,
     double max_r = HALF_SQRT_MAX_DOUBLE, // half for numerical margin error
     double z = std::numeric_limits<double>::min())
   : k_(k), max_r2_( max_r * max_r ), z_(z) {}


  /**
      returns the force constant for pulling
  */
  double get_k() const { return k_; }

  /**
      evaluates the derivative accoridng to the particle position

      @param d core::XYZR position of the particle
      @return a pair with the score and the derivative vector
  */
  std::pair<double, algebra::Vector3D>  eval_deriv(const core::XYZR &d) const {
    double score;
    algebra::Vector3D v_deriv;
    double r2 =  std::pow(d.get_x(), 2) + std::pow(d.get_y(), 2);
    if (r2 > max_r2_ && max_r2_ > 0.0) {
      score = 0.0; 
      v_deriv = algebra::Vector3D(0, 0, 0); // nothing to add to derivatives( = 0);
      return std::make_pair(score, v_deriv);
    }

    double dz = d.get_z() - z_;
    if (dz  > 1e-6) {
      score = k_ * dz;
      v_deriv = algebra::Vector3D(0, 0, k_);
    }
    else if (dz < -1e-6) {
      score = k_ * (-dz);
      v_deriv = algebra::Vector3D(0, 0, -k_);
    }
    else{
      score = 0;
      v_deriv = algebra::Vector3D(0, 0, 0);
    }

    return std::make_pair(score, v_deriv);
  }

  virtual double evaluate_index(Model *m, ParticleIndex pi,
                                DerivativeAccumulator *da) const override
  {
    core::XYZR d(m,pi);
    auto score_deriv = eval_deriv(d);
    auto score = score_deriv.first;
    auto v_deriv = score_deriv.second;
    if (da) {
      IMP_LOG(VERBOSE, "result in " << score
              << " and " << v_deriv << std::endl);
      d.add_to_derivatives(v_deriv, *da);
    }
    return score;
  }

  virtual ModelObjectsTemp do_get_inputs(Model *m,
                                         const ParticleIndexes &pis) const
      override
  { return IMP::get_particles(m, pis);  }

  IMP_SINGLETON_SCORE_METHODS(ZBiasSingletonScore);
  IMP_OBJECT_METHODS(ZBiasSingletonScore);
};

IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_Z_BIAS_SINGLETON_SCORE_H */
