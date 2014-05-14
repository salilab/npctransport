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
#include <IMP/display/Color.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

struct TAMD_chain {
public:
  IMP::Particle* p;
  IMP::Particles centroids;
  IMP::Particles images;
  IMP::Restraints R;
public:
TAMD_chain(IMP::Particle* pp,
           IMP::Particles ccentroids,
           IMP::Particles iimages,
           IMP::Restraints RR):
  p(pp), centroids(ccentroids), images(iimages), R(RR)
  {}
};

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
  IMP::Particle* create_tamd_fg_chain
( SimulationData *sd,
  const ::npctransport_proto::Assignment_FGAssignment &fg_data,
  display::Color c );


/**
  Create a TAMD hierarchy of nlevels depth with a core for
  each d centroids in a lower level, with a real centroid and
  restrained centroid realization

  @param m       Model
  @param pf      A factory for producing singleton particles
                 (the leaves of the chain)
  @param nlevels Number of tamd levels in the hierarchy. If 0 then return a
                 singleton particle.
  @param d       The out degree of each non-leaf node (# of children)
  @param T_factors A list of length nlevels with temeprature at each level
                 from top to bottom
  @param G_factors A list of length nlevels with friction factor (G for gamma)
                 at each level from top to bottom
  @param Ks      Spring constants at each level between a particle and its TAMD
                 image (a list of length nlevels)

  @return a tuple with <root particle, centroids, images, restraints>
          // TODO: explain some more
*/
//boost::tuple<IMP::Particle*, IMP::Particles, IMP::Particles, IMP::Restraints>
TAMD_chain
create_tamd_chain( Model* m,
                   ParticleFactory pf,
                   unsigned int nlevels,
                   unsigned int d,
                   std::vector<double> T_factors,
                    std::vector<double> F_factors,
                   std::vector<double> Ks );

// def _get_ordered_leaves(self, h):



IMPNPCTRANSPORT_END_NAMESPACE

#endif /* IMPNPCTRANSPORT_CREATING_TAMD_PARTICLES_H */
