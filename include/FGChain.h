/**
 * \file creating_tamd_particles.h
 * \brief creating TAMD chains
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_CREATING_TAMD_PARTICLES_H
#define IMPNPCTRANSPORT_CREATING_TAMD_PARTICLES_H

#include "npctransport_config.h"
#include "SimulationData.h"
#include <IMP/base/Object.h>
#include <IMP/display/Color.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

/**
   a chain with a root

    root - root of hierarchy
    beads - fine chain particles
*/
struct Chain : IMP::base::Object {
public:
  IMP::Particle* root;
  IMP::Particles beads;

 public:
  Chain(IMP::Particle* rroot,
        IMP::Particles bbeads,
        std::string name = "chain %1%")
    : base::Object(name),
    root(rroot),
    beads(bbeads)
    { }
};

/** a hierarchy of centroids to represent a polymer chain with
    centroids in the non-leaf nodes, singleton particles in the
    leaves, and TAMD images attached to the centroids by spring

    root - root of hierarchy
    beads - fine chain particles
    centroids - the list of centroids in tree(p)
    images - corresponding TAMD images for each centroid in centroids
    R - corresponding springs that attch each TAMD image to each centroid
*/
struct TAMDChain : public chain{
  typedef chain P;
 public:
  IMP::Particles centroids;
  IMP::Particles images;
  IMP::Restraints R;
 public:
 TAMDChain(IMP::Particle* rroot,
            IMP::Particles fbeads,
            IMP::Particles ccentroids,
            IMP::Particles iimages,
            IMP::Restraints RR,
            std::string name="tamd_chain %1%")
   : P(rroot, fbeads, name),
    centroids(ccentroids),
    images(iimages),
    R(RR)
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

   @param[in,out] sd the simulation data whose model is associated with the
                  new chain. A chain restraint is added to the simulation data
                  scoring object, and the particle is added to the simulation data
                  diffusers list.
   @param[in] fg_data data about the FG chain
   @param[in] c        color of chain particles

   @return chain structure (with root and chain beads)

 */
Chain* create_tamd_fg_chain
( SimulationData *sd,
  const ::npctransport_proto::Assignment_FGAssignment &fg_data,
  display::Color c );


/**
   gets a chain structure from a root of an FG nup
   (by adding its ordered leaves)
*/
Chain* get_chain(atom::Hierarchy root);


IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_CREATING_TAMD_PARTICLES_H */
