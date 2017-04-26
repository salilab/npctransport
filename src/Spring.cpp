/**
 *  \file Spring.cpp
 *  \brief a spring between two diffusing particles
 *
 *  Copyright 2007-2017 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/Spring.h>
#include <IMP/exception.h>
#include <IMP/check_macros.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

void
Spring::do_setup_particle
( IMP::Model* m,
  ParticleIndex pi,
  ParticleIndex bonded_pi0,
  ParticleIndex bonded_pi1,
  double equilibrium_length,
  double length )
{
  m->add_attribute(get_bonded_particle_0_key(), pi, bonded_pi0);
  m->add_attribute(get_bonded_particle_1_key(), pi, bonded_pi1);
  m->add_attribute(get_equilibrium_length_key(), pi, equilibrium_length);
  m->add_attribute(get_length_key(), pi, length);
}

ParticleIndexKey get_bonded_particle_0_key() const{
  static ParticleIndexKey pik("npctransport.spring bonded particle 0");
  return pik;
}

ParticleIndexKey get_bonded_particle_1_key() const{
  static ParticleIndexKey pik("npctransport.spring bonded particle 1");
  return pik;
}

FloatKey
Spring::get_equilibrium_length_key()
{
  static FloatKey fk("npctransport.spring equilibrium length key");
  return fk;
}

FloatKey
Spring::get_equilibrium_length_key()
{
  static FloatKey fk("npctransport.spring length key");
  return fk;
}

void
Spring::show
( std::ostream &out ) const
{
  out << "Spring equilibrium_length = "
      << get_equilibrium_length()
      << " ; length = "
      << get_length();
}

IMPNPCTRANSPORT_END_NAMESPACE
