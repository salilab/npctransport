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
#include <IMP/base/Pointer.h>
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
  Restraints ret = tamd_restraints_;
  ret += FGChain::get_chain_restraints();
  return ret;
}


/************************* declarations ****************/
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


/************* utility methods implementation ***************/

TAMDChain*
create_tamd_chain( ParticleFactory* pf,
                   unsigned int n,
                   unsigned int d,
                   std::vector<double> T_factors,
                    std::vector<double> F_factors,
                   std::vector<double> Ks )
{
  base::Pointer<TAMDChain> ret_chain = new TAMDChain();
  Particles children; // direct children of uppermost centroid
  Model* m=pf->get_model();
  int nlevels = ceil(log2(n));

  if (n==1)
    {
      // Create and return a singleton leaf:
      ret_chain->set_root( pf->create("leaf %1%") );
      return ret_chain.release();
    }

  // Build <d> children recursively and accumulate in ret_chain:
  {
    // factors and n particles in each child chain
    std::vector<double> T_factors1(T_factors.begin()+1, T_factors.end());
    std::vector<double> F_factors1(F_factors.begin()+1, F_factors.end());
    std::vector<double> Ks1(Ks.begin()+1, Ks.end());
    int n1_base = n / d; // per child
    int n_excess = n % d;
    int n_left = n;
    while(n_left > 0){
      int n1 = n1_base + (n_excess-- > 0 ? 1 : 0);
      base::Pointer<TAMDChain> child =
        create_tamd_chain(pf, n1, d, T_factors, F_factors1, Ks1);
      n_left -= n1;
      // Accumulate all results
      children.push_back( child->get_root() );
      ret_chain->centroids += child->centroids;
      ret_chain->images += child->images;
      ret_chain->tamd_restraints_ += child->tamd_restraints_;
      std::cout << "Created chain with " << n1 << " beads, "
                << n_left << " left" << std::endl;
    }
  }
  for(unsigned int i=0 ; i < ret_chain->get_number_of_beads(); i++){
        std::cout << "n = " << n << " Bead # " << i << " - "
                  << ret_chain->get_bead(i) << std::endl;
  }

  // Build centroid of <d> children and store in ret_chain
  {
    std::ostringstream oss;  oss << "centroid " << nlevels;
    Particle* pc =  new Particle( m, oss.str() );
    core::Typed::setup_particle(pc, pf->type_);
    core::XYZR pc_xyzr = core::XYZR::setup_particle(pc);
    pc_xyzr.set_radius(2); // TODO: something smarter?
    pc_xyzr.set_coordinates_are_optimized(true); // TODO: is needed or dangerous
                                                 // - for BD to evaluate it?
    atom::Diffusion pc_Diffusion = atom::Diffusion::setup_particle(pc);
    atom::Hierarchy pc_h = atom::Hierarchy::setup_particle(pc);
    for(unsigned int i = 0; i < children.size(); i++) {
      IMP_LOG_PROGRESS("Adding " << i << "children to pi" <<
                       pc->get_index() << std::endl);
      pc_h.add_child( atom::Hierarchy(children[i]) );
    }
    // base::Pointer<core::ChildrenRefiner> refiner =
    //  new core::ChildrenRefiner( atom::Hierarchy::get_default_traits() );
    //    atom::CenterOfMass::setup_particle(pc, refiner.get());
    atom::CenterOfMass::setup_particle(pc, children);
    ret_chain->set_root( pc );
    ret_chain->centroids.push_back( pc );
  }

  // Build TAMD image of centroid + spring restraint
  {
    std::ostringstream oss;  oss << "TAMD " << nlevels << ".%1%";
    Particle* image = create_tamd_image(  ret_chain->get_root(),
                                          oss.str(),
                                          T_factors[0],
                                          F_factors[0]);
    ret_chain->images.push_back(image);
    base::Pointer<core::HarmonicDistancePairScore> spring=
      new core::HarmonicDistancePairScore(0, Ks[0]);
    ret_chain->tamd_restraints_.push_back
      ( new core::PairRestraint( spring,
                                 ParticlePair(ret_chain->get_root(), image),
                                 oss.str() ) );
    ret_chain->tamd_springs_.push_back( spring );
  }

  return ret_chain.release();
}


/********************** internal utility functions U**************/

/**
   Create a TAMD image of centroid particle p

   @param p_ref reference particle to be tied by spring
   @param name particle name
   @param T_factor temeprature factor
   @param F_factor friction factor

   @return TAMD image particle */
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
  core::Typed::setup_particle( p_ret, core::Typed(p_ref).get_type() );
  IMP::atom::TAMDParticle::setup_particle( p_ret, p_ref, T_factor, F_factor);
  return p_ret;
}




IMPNPCTRANSPORT_END_INTERNAL_NAMESPACE
