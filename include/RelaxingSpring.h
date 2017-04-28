/**
 *  \file IMP/npctransport/RelaxingSpring.h
 *  \brief a spring between two diffusing particles whose resting length
 *         is relaxing towards some equilibrium value
 *
 *  Copyright 2007-2017 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPNPCTRANSPORT_RELAXING_SPRING_H
#define IMPNPCTRANSPORT_RELAXING_SPRING_H

#include "npctransport_config.h"
#include <IMP/Decorator.h>
#include <IMP/decorator_macros.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

//! A decorator for a spring particle connecting two diffusing particles
/** \ingroup helper
    \ingroup decorators
 */

class IMPNPCTRANSPORTEXPORT
RelaxingSpring : public Decorator()
{
  /** Decorate a spring particle that connects two particles
      in no particular order

      @param m the model
      @param pi the particle index
      @param bonded_pi0 particle on first side of spring
      @param bonded_pi1 particle on second side of spring
      @param equilibrium_rest_length the rest length of the spring at equilibrium
      @param rest_length the instantaneous rest length of the spring
  */
  static void do_setup_particle(IMP::Model* m,
                                ParticleIndex pi,
                                ParticleIndex bonded_pi0,
                                ParticleIndex bonded_pi1,
                                double equilibrium_rest_length,
                                double rest_length);


 public:
  IMP_DECORATOR_METHODS(RelaxingSpring, Decorator);


  /** Decorate a spring particke connecting two particles

      @param m the model
      @param pi the particle index
      @param bonded_pi0 particle on first side of spring
      @param bonded_pi1 particle on second side of spring
      @param equilibrium_rest_length the rest length of the spring at equilibrium
      @param rest_length the instantaneous rest length of the spring
  */
  IMP_DECORATOR_SETUP_4(RelaxingSpring,
                        ParticleIndex bonded_pi0,
                        ParticleIndex bonded_pi1,
                        double, equilibrium_rest_length,
                        double, rest_length);

  //! Return true if the particle is an instance of an Transporting
  static bool get_is_setup(Model *m, ParticleIndex pi) {
    return
      m->get_has_attribute(get_equilibrium_rest_length_key(), pi)
      && m->get_has_attribute(get_rest_length_key(), pi);
  }

  ParticleIndexKey get_bonded_particle_0_key() const;

  ParticleIndexKey get_bonded_particle_1_key() const;

  //! get decorator key for spring equilibrium rest length
  static FloatKey get_equilibrium_rest_length_key() const;

  //! get decorator key for spring rest length
  static FloatKey get_rest_length_key() const;

  //! get decorator key for diffusion coefficient of rest length
  static FloatKey get_rest_length_diffusion_coefficient_key() const;

  ParticleIndex get_bonded_particle_0() const{
    return (get_particle()->get_value(get_bonded_particle_0_key()));
  }

  ParticleIndex get_bonded_particle_1() const{
    return (get_particle()->get_value(get_bonded_particle_1_key()));
  }

  IMP_DECORATOR_GET_SET(equilibrium_rest_length, 
			get_equilibrium_rest_length_key(), 
			Float, Float);

  IMP_DECORATOR_GET_SET(rest_length, 
			get_rest_length_key(), 
			Float, Float);

  IMP_DECORATOR_GET_SET(rest_length_diffusion_coefficient,
			get_rest_length_diffusion_coefficient_key(),
			Float, Float);

  void add_to_rest_length_derivative
    ( double d, DerivativeAccumulator* da ) 
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
