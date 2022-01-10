/**
 *  \file IMP/atom/BrownianDynamicsTAMDWithSlabSupport.cpp
 *  \brief Simple molecular dynamics optimizer.
 *
 *  Copyright 2007-2022 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/BrownianDynamicsTAMDWithSlabSupport.h>
#include <IMP/npctransport/RelaxingSpring.h>
#include <IMP/atom/BrownianDynamicsTAMD.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

void
BrownianDynamicsTAMDWithSlabSupport
::do_advance_chunk
(double dtfs,
 double ikt,
 const ParticleIndexes &ps,
 unsigned int begin,
 unsigned int end)
{
  BrownianDynamicsTAMD::do_advance_chunk(dtfs, ikt, ps, begin, end);
  // TODO: update slab here
  // TODO: move to .cpp

  // Relax all FG springs by going over all FGs and then updating all their springs by random diffusion + gradient just as BD of XYZ particles
  // Note that this is inefficient if there are no harmonic springs (backward support), but we only care about performance of new version
  // where they are persent
  double dtfs_ikt(dtfs*ikt);
  for(unsigned int i=0; i<ps.size(); i++){
    if(!RelaxingSpring::get_is_setup(get_model(),
                                     ps[i])) {
      continue;
    }
    RelaxingSpring rs(get_model(), ps[i]);
    double rest_length_derivative(rs.get_rest_length_derivative());
    double rest_length_diffusion_coefficient(rs.get_rest_length_diffusion_coefficient());
    double rest_length(rs.get_rest_length());
    double sigma(std::sqrt(2*dtfs*rest_length_diffusion_coefficient)); // note 2 and not 6 since a single d.o.f.
    rest_length += get_sample(sigma);
    rest_length -= rest_length_derivative*rest_length_diffusion_coefficient*dtfs_ikt;
    rs.set_rest_length(rest_length);
  }
}

IMPNPCTRANSPORT_END_NAMESPACE
