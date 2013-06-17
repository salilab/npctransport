/**
 *  \file SeparateSingletonModifier.cpp
 *  \brief XXXX.
 *
 *  Copyright 2007-8 Sali Lab. All rights reserved.
 *
 */

#include "IMP/npctransport/ExcludeZRangeSingletonScore.h"
#include "IMP/core/XYZR.h"

IMPNPCTRANSPORT_BEGIN_NAMESPACE

ExcludeZRangeSingletonScore::ExcludeZRangeSingletonScore
(double bottom, double top, double k)
: bottom_(bottom),
    top_(top),
    k_(k)
{   }


double ExcludeZRangeSingletonScore::evaluate_index
(Model *m, ParticleIndex pi,
 DerivativeAccumulator *da) const {
  using algebra::Vector3D;
  IMP_OBJECT_LOG;
  core::XYZR d(m, pi);
  if (!d.get_coordinates_are_optimized()) return false;
  // check for violation
  double z = d.get_z();
  double r = d.get_radius();
  double top_violation = top_ - (z - r);
  double bottom_violation = (z + r) - bottom_;
  if(top_violation < 0 || bottom_violation < 0) // out of z-range
    return 0;
  double score= k_* std::min(top_violation, bottom_violation);
  if (da) {
    Vector3D dc(0, 0, -k_); // go to top
    if(bottom_violation < top_violation) // closer to bottom
      dc = -dc;  // go to bottom
    IMP_LOG(VERBOSE, "result in " << score << " and " << dc
            << std::endl);
    d.add_to_derivatives(dc, *da);
  }
  return score;
}

ModelObjectsTemp
ExcludeZRangeSingletonScore::do_get_inputs(Model *m,
                                           const ParticleIndexes &pis) const {
  return IMP::kernel::get_particles(m, pis);
}

IMPNPCTRANSPORT_END_NAMESPACE
