/**
 * \file internal/creating_tamd_particles.h
 * \brief creating TAMD chains internal methods
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 */

#ifndef IMPNPCTRANSPORT_INTERNAL_TAMD_CHAIN_H
#define IMPNPCTRANSPORT_INTERNAL_TAMD_CHAIN_H

#include "../npctransport_config.h"
#include <IMP/npctransport/ParticleFactory.h>
#include <IMP/npctransport/FGChain.h>
#include <IMP/atom/Hierarchy.h>
#include <string>
#include <vector>

IMPNPCTRANSPORT_BEGIN_INTERNAL_NAMESPACE

/** a hierarchy of centroids to represent a polymer chain with
    centroids in the non-leaf nodes, singleton particles in the
    leaves, and TAMD images attached to the centroids by spring

    root - root of hierarchy
    beads - fine chain particles
    centroids - the list of centroids in tree(p)
    images - corresponding TAMD images for each centroid in centroids
    R - corresponding springs that attch each TAMD image to each centroid
*/
struct TAMDChain : public npctransport::FGChain{
  typedef npctransport::FGChain P;
 public:
  IMP::Particles centroids;
  IMP::Particles images;
  IMP::Restraints R;
 public:
  TAMDChain(std::string name = "tamd_chain %1%")
    : P(name)
    {}

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
   Create a TAMD hierarchy of nlevels depth with a core for
   each d centroids in a lower level, with a real centroid and
   restrained centroid realization

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
*/
TAMDChain*
create_tamd_chain( ParticleFactory pf,
                   unsigned int nlevels,
                   unsigned int d,
                   std::vector<double> T_factors,
                   std::vector<double> F_factors,
                   std::vector<double> Ks );

/**
   Create a TAMD image of centroid particle p

   @param p_ref reference particle to be tied by spring
   @param name particle name
   @param T_factor temeprature factor
   @param F_factor friction factor

   @return TAMD image particle
*/
Particle* create_tamd_image( Particle* p_ref,
                             std::string name,
                             double T_factor,
                             double F_factor);

IMPNPCTRANSPORT_END_INTERNAL_NAMESPACE


#endif /* IMPNPCTRANSPORT_INTERNAL_TAMD_CHAIN_H */
