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


double ExcludeZRangeSingletonScore::evaluate
(Particle *p,
 DerivativeAccumulator *da) const {
  using algebra::Vector3D;
  IMP_OBJECT_LOG;
  core::XYZR d(p);
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

ParticlesTemp
ExcludeZRangeSingletonScore::get_input_particles(Particle*p) const {
  return ParticlesTemp(1,p);
}

ContainersTemp
ExcludeZRangeSingletonScore::get_input_containers(Particle*) const {
  return ContainersTemp();
}


void ExcludeZRangeSingletonScore::do_show(std::ostream &out) const {
  out << "ExcludeZRangeSingletonScore " << bottom_ << " to " << top_
      << " ; k = " << k_ << std::endl;
}

IMPNPCTRANSPORT_END_NAMESPACE
