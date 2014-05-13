/**
 * \file creating_particles.h
 * \brief description
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_PARTICLE_FACTORY_H
#define IMPNPCTRANSPORT_PARTICLE_FACTORY_H

#include "npctransport_config.h"
#include "SimulationData.h"
#include <IMP/base/WeakPointer.h>
#include <IMP/core/Typed.h>
#include <IMP/display/Color.h>
#include <string>

IMPNPCTRANSPORT_BEGIN_NAMESPACE


class ParticleFactory {
 public:
  /** The simulation data whose model is associated with new particle
  */
  IMP::base::WeakPointer<SimulationData> sd_;

  /** particle radius (A) */
  double radius_;

  /** diffusion factor */
  double D_factor_;

  /** angular diffusion factor */
  double angular_D_factor_;

  /** particle color */
  display::Color color_;

  /** particle type */
  core::ParticleType type_;

 public:
  /**
     construct a factory that produces particles with specified attributes

    @param sd the simulation data whose model is associated with new particles
    - particles are also saved to the sd diffusers list
    @param radius particle radius (A)
    @param D_factor diffusion factor (relative to that auto-calculated
                    from radius)
    @param angular_D_factor angular diffusion factor (relative to that
                            auto-calculated from radius*D_factor)
    @param color color for new particles
    @param type the type of new particles
   */
 ParticleFactory(SimulationData* sd,
                 double radius,
                 double D_factor,
                 double angular_D_factor,
                 display::Color color,
                 core::ParticleType type)
   : sd_(sd),
    radius_(radius),
    D_factor_(D_factor),
    angular_D_factor_(angular_D_factor),
    color_(color),
    type_(type)
    {}

  /**
     create a particle as specified during construction.
     The particle is decorated with XYZR, RigidBodyDiffusion,
     Color and Typed decorator. In addition, a simulation data
     key is added to it pointing to this->sd_

     @param name - name of particle. If empty, the return value of
     this->type_.get_string() is used
  */
  IMP::Particle* create(std::string name="");

  Model* get_model() {
    return sd_->get_model();
  }
};


IMPNPCTRANSPORT_END_NAMESPACE

#endif  /* IMPNPCTRANSPORT_PARTICLE_FACTORY_H */
