/**
 *  \file SlabWithCylindricalPore.cpp
 *  \brief Decoratr for slab particle with a cylindrical pore
 *
 *  Copyright 2007-2021 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/SlabWithCylindricalPore.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

void
SlabWithCylindricalPore::do_setup_particle(IMP::Model* m,
					   ParticleIndex pi,
					   double thickness,
					   double radius)
{
  SlabWithPore::setup_particle(m, pi, thickness, radius);
  m->add_attribute(get_cylindrical_pore_key(), pi, true);
}

IntKey SlabWithCylindricalPore::get_cylindrical_pore_key(){
  static IntKey k("cylindrical_pore");
  return k;
}

void SlabWithCylindricalPore::show(std::ostream &out) const {
  out << "SlabWithCylindricalPore";
  SlabWithPore::show(out);
}

IMPNPCTRANSPORT_END_NAMESPACE
