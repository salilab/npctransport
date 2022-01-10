/**
 *  \file IMP/npctransport/RelaxingSpring.h
 *  \brief a spring between two diffusing particles whose resting length
 *         is relaxing towards some equilibrium value
 *
 *  Copyright 2007-2022 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPNPCTRANSPORT_RELAXING_SPRING_H
#define IMPNPCTRANSPORT_RELAXING_SPRING_H

#include "npctransport_config.h"
#include <IMP/Decorator.h>
#include <IMP/decorator_macros.h>
#include <IMP/Particle.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

//! A decorator for a spring particle connecting two diffusing particles
/** \ingroup helper
    \ingroup decorators
 */

class IMPNPCTRANSPORTEXPORT
RelaxingSpring : public Decorator
{
  /** Decorate a spring particle that connects two particles
      in no particular order, with a dynamic rest length
      that may relax towards some equilibrium value

      @param m the model
      @param pi the particle index
      @param bonded_pi0 particle on first side of spring
      @param bonded_pi1 particle on second side of spring
      @param equilibrium_rest_length_factor the rest length factor of the spring at equilibrium (relative to sum of the radii of the bonded particles)
      @param rest_length_diffusion_coefficient the diffusion coefficient for the rest length
  */
  static void do_setup_particle(IMP::Model* m,
                                ParticleIndex pi,
                                ParticleIndex bonded_pi0,
                                ParticleIndex bonded_pi1,
                                double equilibrium_rest_length_factor,
                                double rest_length_diffusion_coefficient);


 public:
  IMP_DECORATOR_METHODS(RelaxingSpring, Decorator);


  /** Decorate a spring particle that connects two particles
      in no particular order, with a dynamic rest length
      that may relax towards some equilibrium value

      @param m the model
      @param pi the particle index
      @param bonded_pi0 particle on first side of spring
      @param bonded_pi1 particle on second side of spring
      @param equilibrium_rest_length_factor the rest length factor of the spring at equilibrium (relative to sum of the radii of the bonded particles)
      @param rest_length_diffusion_coefficient the diffusion
             coefficient for the rest length
       */
  IMP_DECORATOR_SETUP_4(RelaxingSpring,
                        ParticleIndex, bonded_pi0,
                        ParticleIndex, bonded_pi1,
                        double, equilibrium_rest_length_factor,
                        double, rest_length_diffusion_coefficient);

  //! Return true if the particle is an instance of an Transporting
  static bool get_is_setup(Model *m, ParticleIndex pi) {
    return
      m->get_has_attribute(get_equilibrium_rest_length_factor_key(), pi)
      && m->get_has_attribute(get_rest_length_key(), pi);
  }

  static ParticleIndexKey get_bonded_particle_0_key();

  static ParticleIndexKey get_bonded_particle_1_key();

  //! get decorator key for spring equilibrium rest length factor
  static FloatKey get_equilibrium_rest_length_factor_key();

  //! get decorator key for spring rest length
  static FloatKey get_rest_length_key();

  //! get decorator key for diffusion coefficient of rest length
  static FloatKey get_rest_length_diffusion_coefficient_key();

  Particle* get_bonded_particle_0() const{
    Particle* this_p= get_particle();
    return this_p->get_value(get_bonded_particle_0_key());
  }

  Particle* get_bonded_particle_1() const{
    Particle* this_p= get_particle();
    return this_p->get_value(get_bonded_particle_1_key());
  }

  ParticleIndex get_bonded_particle_index_0() const{
    return get_bonded_particle_0()->get_index();
  }

  ParticleIndex get_bonded_particle_index_1() const{
    return get_bonded_particle_1()->get_index();
  }


  IMP_DECORATOR_GET_SET(equilibrium_rest_length_factor,
			get_equilibrium_rest_length_factor_key(),
			Float, Float);

  IMP_DECORATOR_GET_SET(rest_length,
			get_rest_length_key(),
			Float, Float);

  IMP_DECORATOR_GET_SET(rest_length_diffusion_coefficient,
			get_rest_length_diffusion_coefficient_key(),
			Float, Float);

  void add_to_rest_length_derivative
    ( double d, DerivativeAccumulator& da )
  {
    get_model()->add_to_derivative( get_rest_length_key(),
				    get_particle_index(),
				    d,
				    da );
  }

  double get_rest_length_derivative()
  {
    return get_model()->get_derivative( get_rest_length_key(),
					get_particle_index() );
  }

};



IMP_DECORATORS(RelaxingSpring, RelaxingSprings, IMP::Decorators);

IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_RELAXING_SPRING_H */
