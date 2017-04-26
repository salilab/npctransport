/**
 *  \file IMP/npctransport/Spring.h
 *  \brief A decorator for a spring particle connecting two diffusing particles
 *
 *  Copyright 2007-2017 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPNPCTRANSPORT_SPRING_H
#define IMPNPCTRANSPORT_SPRING_H

#include "npctransport_config.h"
#include <IMP/Decorator.h>
#include <IMP/decorator_macros.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

//! A decorator for a spring particle connecting two diffusing particles
/** \ingroup helper
    \ingroup decorators
 */

class IMPNPCTRANSPORTEXPORT
Spring : public Decorator()
{
  /** Decorate a spring particle that connects two particles
      in no particular order

      @param m the model
      @param pi the particle index
      @param bonded_pi0 particle on first side of spring
      @param bonded_pi1 particle on second side of spring
      @param equilibrium_length the length of the spring at equilibrium
      @param length the instantaneous length of the spring
  */
  static void do_setup_particle(IMP::Model* m,
                                ParticleIndex pi,
                                ParticleIndex bonded_pi0,
                                ParticleIndex bonded_pi1,
                                double equilibrium_length,
                                double length);


 public:
  IMP_DECORATOR_METHODS(Spring, Decorator);


  /** Decorate a spring particke connecting two particles

      @param m the model
      @param pi the particle index
      @param bonded_pi0 particle on first side of spring
      @param bonded_pi1 particle on second side of spring
      @param equilibrium_length the length of the spring at equilibrium
      @param length the instantaneous length of the spring
  */
  IMP_DECORATOR_SETUP_4(Spring,
                        ParticleIndex bonded_pi0,
                        ParticleIndex bonded_pi1,
                        double, equilibrium_length,
                        double, length);

  //! Return true if the particle is an instance of an Transporting
  static bool get_is_setup(Model *m, ParticleIndex pi) {
    return
      m->get_has_attribute(get_equilibrium_length_key(), pi)
      && m->get_has_attribute(get_length_key(), pi);
  }

  ParticleIndex get_bonded_particle_0() const{
    return (get_particle()->get_value(get_bonded_particle_0_key()));
  }

  ParticleIndexKey get_bonded_particle_0_key() const;

  ParticleIndex get_bonded_particle_1() const{
    return (get_particle()->get_value(get_bonded_particle_1_key()));
  }

  ParticleIndexKey get_bonded_particle_1_key() const;

  //! sets the equilibrium length of the spring
  void set_equilibrium_length(double equilibrium_length) {
    get_particle()->set_value(get_equilibrium_length_key(),
                              equilibrium_length);
  }

  //! returns the equilibrium length of the spring
  double get_equilibrium_length() const {
    return (get_particle()->get_value(get_equilibrium_length_key()));
  }

  //! get decorator key for spring equilibrium length
  static FloatKey get_equilibrium_length_key() const;

  //! sets the length of the spring
  void set_length(double length) {
    get_particle()->set_value(get_length_key(),
                              length);
  }

  //! returns the length of the spring
  double get_length() const {
    return (get_particle()->get_value(get_length_key()));
  }

  //! get decorator key for spring length
  static FloatKey get_length_key() const;

};



IMP_DECORATORS(Spring, Springs, IMP::Decorators);

IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_SPRING_H */
