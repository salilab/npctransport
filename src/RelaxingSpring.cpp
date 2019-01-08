/**
 *  \file RelaxingSpring.cpp
 *  \brief a spring between two diffusing particles whose resting length
 *         is relaxing towards some equilibrium value
 *
 *  Copyright 2007-2019 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/RelaxingSpring.h>
#include <IMP/exception.h>
#include <IMP/check_macros.h>
#include <IMP/core/XYZR.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

void
RelaxingSpring::do_setup_particle
( IMP::Model* m,
  ParticleIndex pi,
  ParticleIndex bonded_pi0,
  ParticleIndex bonded_pi1,
  double equilibrium_rest_length_factor,
  double rest_length_diffusion_coefficient )
{
  IMP_ALWAYS_CHECK(IMP::core::XYZR::get_is_setup(m, bonded_pi0) &&
                   IMP::core::XYZR::get_is_setup(m, bonded_pi1),
                   "can only bond particles that are of xyzr type",
                   IMP::ValueException);
  m->add_attribute(get_bonded_particle_0_key(), pi,
		   bonded_pi0);
  m->add_attribute(get_bonded_particle_1_key(), pi,
		   bonded_pi1);
  m->add_attribute(get_equilibrium_rest_length_factor_key(), pi,
		   equilibrium_rest_length_factor);
  IMP::core::XYZR xyzr0(m, bonded_pi0);
  IMP::core::XYZR xyzr1(m, bonded_pi1);
  double r0= xyzr0.get_radius();
  double r1= xyzr0.get_radius();
  double init_rest_length= equilibrium_rest_length_factor * (r0+r1);
  m->add_attribute(get_rest_length_key(), pi,
		   init_rest_length); // initially at equilibrium by default
  m->add_attribute(get_rest_length_diffusion_coefficient_key(), pi,
		   rest_length_diffusion_coefficient);
}

ParticleIndexKey
RelaxingSpring::get_bonded_particle_0_key(){
  static ParticleIndexKey pik("npctransport.spring bonded particle 0");
  return pik;
}

ParticleIndexKey
RelaxingSpring::get_bonded_particle_1_key(){
  static ParticleIndexKey pik("npctransport.spring bonded particle 1");
  return pik;
}

FloatKey
RelaxingSpring::get_equilibrium_rest_length_factor_key()
{
  static FloatKey fk("npctransport.spring equilibrium rest length factor key");
  return fk;
}

FloatKey
RelaxingSpring::get_rest_length_key()
{
  static FloatKey fk("npctransport.spring rest_length key");
  return fk;
}

FloatKey
RelaxingSpring::get_rest_length_diffusion_coefficient_key()
{
  static FloatKey fk("npctransport.spring rest length diffusion coefficient key");
  return fk;
}

void
RelaxingSpring::show
( std::ostream &out ) const
{
  out << "RelaxingSpring equilibrium rest length factor = "
      << get_equilibrium_rest_length_factor()
      << "; rest length = "
      << get_rest_length()
      << "; rest length diffusion_coefficient = "
      << get_rest_length_diffusion_coefficient();
}

IMPNPCTRANSPORT_END_NAMESPACE
