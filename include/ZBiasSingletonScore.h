/**
 *  \file ZBiasSingletonScore.h
 *  \brief score that biases particles to go down the Z axis
 *
 *  Copyright 2007-8 Sali Lab. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_Z_BIAS_SINGLETON_SCORE_H
#define IMPNPCTRANSPORT_Z_BIAS_SINGLETON_SCORE_H

#include "npctransport_config.h"
#include <IMP/SingletonScore.h>
#include <IMP/singleton_macros.h>
#include <IMP/algebra/Vector3D.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

class IMPNPCTRANSPORTEXPORT ZBiasSingletonScore
    : public SingletonScore {
 private:
  algebra::Vector3D v_deriv_; // deriv of force, which is of form (0,0,k)
 public:
  /**
     Exclude particles from the range of z coordinates [bottom_..top_]
     with repulsive force constant k
   */
  ZBiasSingletonScore(double k)
    { set_k(k); }


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
  { return IMP::kernel::get_particles(m, pis);  }

  IMP_SINGLETON_SCORE_METHODS(ZBiasSingletonScore);
  IMP_OBJECT_METHODS(ZBiasSingletonScore);
};

IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_Z_BIAS_SINGLETON_SCORE_H */
