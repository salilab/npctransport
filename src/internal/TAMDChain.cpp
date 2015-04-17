/**
 * \file internal/creating_tamd_particles.cpp
 * \brief creating TAMD chains internal methods
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 */

#include <IMP/npctransport/internal/TAMDChain.h>
#include <IMP/atom/Diffusion.h>
#include <IMP/atom/Hierarchy.h>
#include <IMP/atom/CenterOfMass.h>
#include <IMP/atom/Mass.h>
#include <IMP/atom/TAMDParticle.h>
#include <IMP/Pointer.h>
#include <IMP/core/ChildrenRefiner.h>
#include <IMP/core/PairRestraint.h>
#include <IMP/core/rigid_bodies.h>
#include <IMP/core/XYZR.h>
#include <IMP/display/Colored.h>

#include <boost/tuple/tuple.hpp>
#include <cmath>
#include <cstdio>
#include <sstream>

IMPNPCTRANSPORT_BEGIN_INTERNAL_NAMESPACE

/********************* TAMDChain methods ****************/

Restraints TAMDChain::get_chain_restraints()
{
  IMP_USAGE_CHECK( is_initialized,
                   "TAMD chain must be initialized by restraint generation");
  Restraints ret = tamd_restraints_;
  ret += FGChain::get_chain_restraints();
  return ret;
}


// add d children and n beads to root root_h (accessory method for
// create_tamd_chain()) ; params are same as create_tamd_chain() ;
// Set root as centroid of its children (with updated list)
void TAMDChain::add_children
( ParticleFactory* pf,
  unsigned int n,
  unsigned int d,
  std::vector<double> T_factors,
  std::vector<double> F_factors,
  std::vector<double> Ks )
{
  std::vector<double> T_factors1(T_factors.begin()+1, T_factors.end());
  std::vector<double> F_factors1(F_factors.begin()+1, F_factors.end());
  std::vector<double> Ks1(Ks.begin()+1, Ks.end());
  int n1_base = n / d; // baseline for beads per child
  int n_excess = n % d;
  int n_left = n;
  while(n_left > 0){
    int n1 = n1_base + (n_excess-- > 0 ? 1 : 0); // actual beads per child
    Pointer<TAMDChain> child_chain =
      create_tamd_chain(pf, n1, d, T_factors, F_factors1, Ks1);
    n_left -= n1;
    get_root().add_child( child_chain->get_root() );
    centroids_ += child_chain->centroids_;
    images_ += child_chain->images_;
    tamd_restraints_ += child_chain->tamd_restraints_;
    tamd_springs_ += child_chain->tamd_springs_;
    std::cout << "Added child chain with " << n1 << " beads; "
              << n_left << " left" << std::endl;
  }

  // Setup root as center of mass of its children:
  // Pointer<core::ChildrenRefiner> refiner =
  //  new core::ChildrenRefiner( atom::Hierarchy::get_default_traits() );
  atom::CenterOfMass::setup_particle(get_root(),
                                     get_root().get_children());

}


/************************* declarations ****************/


/**
   Create a TAMD image of centroid particle p
   ( accessory function for create_tamd_chain() )

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



/************* factory functions implementation ***************/

TAMDChain*
create_tamd_chain( ParticleFactory* pf,
                   unsigned int n,
                   unsigned int d,
                   std::vector<double> T_factors,
                    std::vector<double> F_factors,
                   std::vector<double> Ks )
{
  // Exceptionalize singletons (recursion stop condition)
  if (n==1)
    {
      return create_singleton_tamd_chain( pf );
    }

  // Initial preparation (can't use log2(n) since that requires C99):
  unsigned int nlevels = ceil(log(n) / log(2.));

  // Setup root centroid with xyzr, diffusion, type and mass
  std::ostringstream root_name_oss;
  root_name_oss << pf->type_.get_string()
                << " centroid " << nlevels << " %1%";
  Particle* root = new Particle
    ( pf->get_model(), root_name_oss.str() );
  atom::Hierarchy root_h = atom::Hierarchy::setup_particle(root);
  core::XYZR root_xyzr = core::XYZR::setup_particle(root);
  // note: radius will affect diffusion coefficient + visual
  //       sqrt(n) to approximate a coil
  root_xyzr.set_radius(pf->get_radius() * sqrt(n) );
  root_xyzr.set_coordinates_are_optimized(true); // TODO: needed or dangerous
                                                   // for BD to evaluate it?
  IMP_LOG_PROGRESS("root radius " << root_xyzr.get_radius() << std::endl);
  atom::Diffusion::setup_particle(root); // TODO: is needed?
  core::Typed::setup_particle(root,
                              core::ParticleType("TAMD Centroid") );
  atom::Mass::setup_particle(root, 1.0); // dummy - will be updated by CenterOfMass

  // Build TAMD image of root + tamd spring restraint:
  std::string image_name = "Image " + root->get_name();
  Particle* image = create_tamd_image(  root,
                                        image_name,
                                        T_factors[0],
                                        F_factors[0]);
  Pointer<core::HarmonicDistancePairScore> tamd_spring=
    new core::HarmonicDistancePairScore(0, Ks[0]);
  Pointer<IMP::Restraint> tamd_restraint = new core::PairRestraint
    ( tamd_spring, ParticlePair(root, image), image_name );

  // build TAMD chain object with children and return it:
  IMP_NEW(TAMDChain, ret_chain, ());
  ret_chain->set_root(root);
  ret_chain->centroids_.push_back( root_h );
  ret_chain->images_.push_back(image);
  ret_chain->tamd_restraints_.push_back(tamd_restraint);
  ret_chain->tamd_springs_.push_back( tamd_spring );
  ret_chain->add_children(pf, n, d, T_factors, F_factors, Ks);
  ret_chain->is_initialized = true;
  for(unsigned int i=0 ; i < ret_chain->get_number_of_beads(); i++){
        std::cout << "n = " << n << " Bead # " << i << " - "
                  << ret_chain->get_bead(i) << std::endl;
  }
  return ret_chain.release();
}


//! create a TAMDChain with a single bead (that is also the root),
//! which is produced by pf (accessory function for create_tamd_chain())
TAMDChain* create_singleton_tamd_chain( ParticleFactory* pf)
{
  IMP_NEW(TAMDChain, ret_chain, ());
  Particle* singleton = pf->create("leaf %1%") ;
  ret_chain->set_root( singleton );
  ret_chain->is_initialized = true;
  return ret_chain.release();
}


/********************** internal accessory functions **************/


/**
   Create a TAMD image of centroid particle p

   @param p_ref reference particle to be tied by spring
   @param name particle name
   @param T_factor temeprature factor
   @param F_factor friction factor

   @return TAMD image particle */
// TODO: make a private method? or add factory class and consolidate
//       all factory stuff there? (eg add_children, etc.)
Particle* create_tamd_image( Particle* p_ref,
                             std::string name,
                             double T_factor,
                             double F_factor){
  //TAMD image of centroid
  Particle* p_ret = new IMP::Particle(p_ref->get_model(), name);
  IMP::core::XYZR p_ret_xyzr = IMP::core::XYZR::setup_particle( p_ret );
  IMP::core::XYZR p_ref_xyzr = IMP::core::XYZR(p_ref);
  p_ret_xyzr.set_coordinates( p_ref_xyzr.get_coordinates() );
  p_ret_xyzr.set_radius( p_ref_xyzr.get_radius() );
  p_ret_xyzr.set_coordinates_are_optimized( true );
  IMP::atom::Hierarchy::setup_particle( p_ret );
  IMP::atom::Diffusion::setup_particle( p_ret ); // diffusion coefficient?!
  IMP::atom::Mass::setup_particle( p_ret, atom::Mass(p_ref).get_mass() );
  std::string type_name = "TAMD Image " +
    core::Typed(p_ref).get_type().get_string();
  core::Typed::setup_particle( p_ret, core::ParticleType(type_name) );
  IMP::atom::TAMDParticle::setup_particle( p_ret, p_ref, T_factor, F_factor);
  return p_ret;
}








IMPNPCTRANSPORT_END_INTERNAL_NAMESPACE
