/**
 *  \file SlabWithCylindricalPore.cpp
 *  \brief Decoratr for slab particle with a cylindrical pore
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/SlabWithPore.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

void 
SlabWithPore::do_setup_particle(IMP::Model* m,
				ParticleIndex pi,
				double thickness, 
				double radius)
{
  m->add_attribute(get_thickness_key(), pi, thickness, false);
  m->add_attribute(get_radius_key(), pi, r, false);
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
      << " ; radius=" << get_radius();
}

IMPNPCTRANSPORT_END_NAMESPACE
