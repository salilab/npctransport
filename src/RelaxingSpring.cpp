/**
 *  \file RelaxingSpring.cpp
 *  \brief a spring between two diffusing particles whose resting length
 *         is relaxing towards some equilibrium value
 *
 *  Copyright 2007-2017 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/RelaxingSpring.h>
#include <IMP/exception.h>
#include <IMP/check_macros.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

void
RelaxingSpring::do_setup_particle
( IMP::Model* m,
  ParticleIndex pi,
  ParticleIndex bonded_pi0,
  ParticleIndex bonded_pi1,
  double equilibrium_rest_length,
  double rest_length_diffusion_coefficient )
{
  m->add_attribute(get_bonded_particle_0_key(), pi, 
		   bonded_pi0);
  m->add_attribute(get_bonded_particle_1_key(), pi, 
		   bonded_pi1);
  m->add_attribute(get_equilibrium_rest_length_key(), pi, 
		   equilibrium_rest_length);
  m->add_attribute(get_rest_length_key(), pi, 
		   equilibrium_rest_length); // initially equal by default
  m->add_attribute(get_rest_length_diffusion_coefficient_key(), pi, 
		   rest_length_diffusion_coefficient);
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
RelaxingSpring::get_equilibrium_rest_length_key() const
{
  static FloatKey fk("npctransport.spring equilibrium rest_length key");
  return fk;
}

FloatKey
RelaxingSpring::get_rest_length_key() const
{
  static FloatKey fk("npctransport.spring rest_length key");
  return fk;
}

FloatKey
RelaxingSpring::get_rest_length_diffusion_coefficient_key() const
{
  static FloatKey fk("npctransport.spring rest length diffusion coefficient key");
  return fg;
}

void
RelaxingSpring::show
( std::ostream &out ) const
{
  out << "RelaxingSpring equilibrium_rest_length = "
      << get_equilibrium_rest_length()
      << " ; rest_length = "
      << get_rest_length();
}

IMPNPCTRANSPORT_END_NAMESPACE
