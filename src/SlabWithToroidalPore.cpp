/**
 *  \file SlabWithToroidalPore.cpp
 *  \brief Decoratr for slab particle with a toroidal pore
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/SlabWithToroidalPore.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

void 
SlabWithToroidalPore::do_setup_particle(IMP::Model* m,
						ParticleIndex pi,
						double thickness, 
						double radius)
{
  SlabWithPore::do_setup_particle(m, pi, thickness, radius);
  m->add_attribute(get_toroidal_pore_key(), pi, true, false);
}

IntKey SlabWithToroidalPore::get_toroidal_pore_key(){
  static IntKey k("toroidal_pore");
  return k;
}

void SlabWithToroidalPore::show(std::ostream &out) const {
  out << "SlabWithToroidalPore";
  SlabWithPore::show(out);
}

IMPNPCTRANSPORT_END_NAMESPACE
