/**
 * \file creating_particles.h
 * \brief description
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_CREATING_PARTICLES_H
#define IMPNPCTRANSPORT_CREATING_PARTICLES_H

#include "npctransport_config.h"
#include "SimulationData.h"
#include <IMP/display/Color.h>
#include <IMP/core/Typed.h>
#include <IMP/container/PairContainerSet.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

/**
   create a particle associated with the model of sd, with the given params

   @param sd the simulation data whose model is associated with the new particle
             - the particle is also saved to the sd diffusers list
   @param radius particle radius (A)
   @param angular_D_factor angular diffusion factor
   @param D_factor diffusion factor
   @param c color for the particle
   @param type the type of the particle
   @param name particle name
 */
IMPNPCTRANSPORTEXPORT Particle* create_particle(SimulationData *sd,
                                                double radius,
                                                double angular_D_factor,
                                                double D_factor,
                                                display::Color c,
                                                core::ParticleType type,
                                                std::string name);

/**
   Create a chain particle hierarchy, associated with the model of sd

   @param sd[out] the simulation data whose model is associated with the
                  new chain. A chain restraint is added to the simulation data
                  itself, and the particle is added to the simulation data
                  diffusers list.
   @param n       number of particles in chain
   @param radius  the radius of each particle
   @param angular_D_factor   angular diffusion factor of chain particles
   @param D_factor diffusion factor of chain particles
   @param rest_length_factor the spring resting distance between consecutive
                             particles is [rest_length_factor * 2 * radius]
   @param spring_constant    the spring constant restraining consecutive chain
                             particles
   @param c        color of chain particles
   @param t        the type of the particles in the chain
   @param name     name of the particle at the root of the chain hierarchy
   @param bonds[out] the list of consecutive bonded particles is added as
                     a container in bonds

   @return chain particle (a particle that is the root of the chain hierarchy)

 */
IMPNPCTRANSPORTEXPORT Particle*
create_chain(SimulationData *sd, int n, double radius,
             double angular_D_factor,
             double D_factor, double rest_length_factor,
             double spring_constant,
             display::Color c,
             core::ParticleType t, std::string name);

IMPNPCTRANSPORT_END_NAMESPACE

#endif  /* IMPNPCTRANSPORT_CREATING_PARTICLES_H */
