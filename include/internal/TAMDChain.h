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
#include <IMP/core/DistancePairScore.h>
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
    springs - springs stored in R (in same order)
*/
class TAMDChain : public npctransport::FGChain{
  typedef npctransport::FGChain P;
 public:
  // TODO: keep those in the hierarchy?
  IMP::Particles centroids;
  IMP::Particles images;
 private:
  IMP::Restraints tamd_restraints_;
  core::HarmonicDistancePairScores tamd_springs_;

 public:
  /** Initialize a TAMD Chain that still doesn't have a root or beads,
      with backbone k and rest length factor as specified, to be utilized
      once an actual root with bead leafs is added
  */
  TAMDChain(
           double backbone_k = 0.0,
           double rest_length_factor = 1.0,
           std::string name = "tamd_chain %1%")
    : P(nullptr, backbone_k, rest_length_factor, name)
    {}

  /**
     initialized a TAMDChain whose root is root, with
     centroids images and tamd_restraints as specified, and backone k
     and rest length factor as specified
   */
 TAMDChain(IMP::Particle* root,
           IMP::Particles ccentroids,
           IMP::Particles iimages,
           IMP::Restraints tamd_restraints,
           IMP::core::HarmonicDistancePairScores tamd_springs,
           double backbone_k = 0.0,
           double rest_length_factor = 1.0,
           std::string name="tamd_chain %1%")
   : P(root, backbone_k, rest_length_factor, name),
    centroids(ccentroids),
    images(iimages),
    tamd_restraints_(tamd_restraints),
    tamd_springs_(tamd_springs)
    { }


  /**
      return the restraint accosciated with internal interactions by this chain.
      In particular for a TAMD chain, this includes both the bonds between beads,
      and the TAMD spring restraints.
  */
  virtual Restraints get_chain_restraints() IMP_OVERRIDE;


  IMP_OBJECT_METHODS(TAMDChain);

 /******************* utility methods ********************/

  friend TAMDChain*
    create_tamd_chain( ParticleFactory* pf,
                       unsigned int nlevels,
                       unsigned int d,
                       std::vector<double> T_factors,
                       std::vector<double> F_factors,
                       std::vector<double> Ks );

};


/************* utility methods ***************/

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
create_tamd_chain( ParticleFactory* pf,
                   unsigned int nlevels,
                   unsigned int d,
                   std::vector<double> T_factors,
                   std::vector<double> F_factors,
                   std::vector<double> Ks );




IMPNPCTRANSPORT_END_INTERNAL_NAMESPACE


#endif /* IMPNPCTRANSPORT_INTERNAL_TAMD_CHAIN_H */
