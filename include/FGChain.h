/**
 * \file FGChain.h
 * \brief creating TAMD chains
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_CREATING_TAMD_PARTICLES_H
#define IMPNPCTRANSPORT_CREATING_TAMD_PARTICLES_H

#include "npctransport_config.h"
#include "SimulationData.h"
#include <IMP/base/Object.h>
#include <IMP/base/nullptr.h>
#include <IMP/display/Color.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

/**
   a chain with a root

    root - root of hierarchy
    beads - fine chain particles
*/
struct FGChain : IMP::base::Object {
public:
  IMP::Particle* root;
  IMP::Particles beads;

 public:
 FGChain(std::string name = "chain %1%")
   : base::Object(name),
    root(nullptr) {}

  FGChain(IMP::Particle* rroot,
        IMP::Particles bbeads,
        std::string name = "chain %1%")
    : base::Object(name),
    root(rroot),
    beads(bbeads)
    { }
};


/**
   Create a chain particle hierarchy, associated with the model of sd,
   with restraint bonding consecutive particles, according to the
   parameters specified in fg_data.

   Notes:

   All nodes in the hierarchy share the same type, based on the type
   string fg_data.type().

   The rest length between two consecutive chain beads is
   fg_data.radius() * 2.0 * fg_data.rest_length_factor() and the
   spring constant is the simulation backbone_k parameter.

   If fg_data.is_tamd() is true, created a TAMD hierarchy, otherwise
   a simple parent + beads structure

   @param[in,out] sd the simulation data whose model is associated with the
                  new chain. A chain restraint is added to the simulation data
                  scoring object, and the particle is added to the simulation data
                  diffusers list.
   @param[in] fg_data data about the FG chain
   @param[in] c        color of chain particles

   @return chain structure (with root and chain beads)

 */
FGChain* create_fg_chain
( SimulationData *sd,
  const ::npctransport_proto::Assignment_FGAssignment &fg_data,
  display::Color c );


/**
   gets a chain structure from a root of an FG nup
   (by adding its ordered leaves)
*/
FGChain* get_fg_chain(atom::Hierarchy root);

/**
   gets a chain structure from a root of an FG nup
   (by adding its ordered leaves)
*/
inline FGChain* get_fg_chain(Particle* p_root);

//// IMPLEMENTATION

FGChain* get_fg_chain(Particle* p_root)
{
  IMP_USAGE_CHECK(atom::Hierarchy::get_is_setup(p_root),
                  "p_root must be of type hierarchy");
  return get_fg_chain(atom::Hierarchy(p_root));
};


IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_CREATING_TAMD_PARTICLES_H */
