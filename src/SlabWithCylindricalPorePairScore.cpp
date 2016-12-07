/**
 *  \file SlabWithCylindricalPorePairScore.cpp
 *  \brief XXXX.
 *
 *  Copyright 2007-8 Sali Lab. All rights reserved.
 *
 */

#include "IMP/npctransport/SlabWithCylindricalPorePairScore.h"

IMPNPCTRANSPORT_BEGIN_NAMESPACE

SlabWithCylindricalPorePairScore::SlabWithCylindricalPorePairScore
(double k)
    : k_(k)
{
    IMP_LOG_PROGRESS("Constructing a slab with cylindrical pore singleton score"
		     << "; k=" << k << std::endl);

}

ModelObjectsTemp SlabWithCylindricalPorePairScore::do_get_inputs(
    Model *m, const ParticleIndexes &pis) const {
  return IMP::get_particles(m, pis);
}

algebra::Vector3D
SlabWithCylindricalPorePairScore::get_displacement_direction
(SlabWithCylindricalPore const& slab, const algebra::Vector3D &v) const
{
  update_cached_slab_params(slab);
  return get_displacement_vector(v).second;
}

double
SlabWithCylindricalPorePairScore::get_displacement_magnitude
(SlabWithCylindricalPore const&slab, const algebra::Vector3D &v) const
{
  update_cached_slab_params(slab);
  return get_displacement_vector(v).first;
}



IMPNPCTRANSPORT_END_NAMESPACE
