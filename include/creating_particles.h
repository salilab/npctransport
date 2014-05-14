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
   Create a chain particle hierarchy, associated with the model of sd,
   with restraint bonding consecutive particles, according to the
   parameters specified in fg_data.

   The rest length between two consecutive chain beads is
   fg_data.radius() * 2.0 * fg_data.rest_length_factor() and the
   spring constant is the simulation backbone_k parameter.

   @param[in,out] sd the simulation data whose model is associated with the
                  new chain. A chain restraint is added to the simulation data
                  scoring object, and the particle is added to the simulation data
                  diffusers list.
   @param[in] fg_data data about the FG chain
   @param[in] c        color of chain particles

   @return chain particle (a particle that is the root of the chain hierarchy)

 */
IMPNPCTRANSPORTEXPORT Particle *create_fg_chain(
    SimulationData *sd,
    const ::npctransport_proto::Assignment_FGAssignment &fg_data,
    display::Color c);

IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_CREATING_PARTICLES_H */
