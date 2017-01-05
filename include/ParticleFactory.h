/**
 * \file ParticleFactory.h
 * \brief description
 *
 *  Copyright 2007-2017 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_PARTICLE_FACTORY_H
#define IMPNPCTRANSPORT_PARTICLE_FACTORY_H

#include "npctransport_config.h"
#include "SimulationData.h"
#include <IMP/WeakPointer.h>
#include <IMP/core/Typed.h>
#include <IMP/display/Color.h>

#include <string>

IMPNPCTRANSPORT_BEGIN_NAMESPACE


class IMPNPCTRANSPORTEXPORT ParticleFactory : public IMP::Object {
 public:
  /** The simulation data whose model is associated with new particle
  */
  IMP::WeakPointer<SimulationData> sd_;

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
     construct a factory that produces diffusing particles with specified attributes,
     with a default mass of 1.0

    @param sd the simulation data whose model is associated with new particles
    - particles are decorated with a simulation data attribute to mark their owner
    @param radius particle radius (A)
    @param D_factor diffusion factor (relative to that auto-calculated
                    from radius). If 0.0, no diffusion or angular diffusion
                    is set up.
    @param angular_D_factor angular diffusion factor (relative to that
                            auto-calculated from radius times D_factor). If
                            non-positive, do not setup angular rigid
                            body diffusion (still set up Diffusion if
                            D_factor>0.0)
    @param color color for new particles
    @param type the type of new particles
    @param name object name
   */
 ParticleFactory(SimulationData* sd,
                 double radius,
                 double D_factor,
                 double angular_D_factor,
                 display::Color color,
                 core::ParticleType type,
                 std::string name = "Particle factory %1%")
   :
  IMP::Object(name),
    sd_(sd),
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

  //! return model associated with this factory
  Model* get_model() {
    return sd_->get_model();
  }

  //! return SimulationData object associated with this factory
  SimulationData* get_simulation_data() {
    return sd_;
  }

  //! return radius of generated particles
  double get_radius() const { return radius_; }


  IMP_OBJECT_METHODS(ParticleFactory);
};


IMPNPCTRANSPORT_END_NAMESPACE

#endif  /* IMPNPCTRANSPORT_PARTICLE_FACTORY_H */
