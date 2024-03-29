/**
 *  \file SlabWithCylindricalPore.cpp
 *  \brief Decoratr for slab particle with a cylindrical pore
 *
 *  Copyright 2007-2022 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/SlabWithPore.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

void
SlabWithPore::do_setup_particle(IMP::Model* m,
				ParticleIndex pi,
				double thickness,
				double pore_radius)
{
  m->add_attribute(get_thickness_key(), pi, thickness, false/*is_optimizable*/);
  m->add_attribute(get_pore_radius_key(), pi, pore_radius, true/*is_optimizable*/);
  SlabWithPore(m,pi).set_pore_radius_is_optimized(false); //default
}

FloatKey SlabWithPore::get_thickness_key() {
  static FloatKey fk("thickness");
  return fk;
}

FloatKey SlabWithPore::get_pore_radius_key() {
  static FloatKey fk("pore_radius");
  return fk;
}

void SlabWithPore::show(std::ostream &out) const {
  out << "SlabWithPore thickness="
      << get_thickness()
      << " ; radius=" << get_pore_radius();
}

IMPNPCTRANSPORT_END_NAMESPACE
