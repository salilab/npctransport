/**
 * \file internal/creating_tamd_particles.cpp
 * \brief creating TAMD chains internal methods
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 */

#include <IMP/npctransport/internal/creating_tamd_particles.h>
#include <IMP/atom/Diffusion.h>
#include <IMP/atom/Hierarchy.h>
#include <IMP/atom/CenterOfMass.h>
#include <IMP/atom/Mass.h>
#include <IMP/atom/TAMDParticle.h>
#include <IMP/base/Pointer.h>
#include <IMP/core/ChildrenRefiner.h>
#include <IMP/core/rigid_bodies.h>
#include <IMP/core/XYZR.h>
#include <IMP/display/Colored.h>

#include <boost/tuple/tuple.hpp>
#include <cmath>
#include <cstdio>
#include <sstream>

IMPNPCTRANSPORT_BEGIN_INTERNAL_NAMESPACE


TAMD_chain*
create_tamd_chain( ParticleFactory pf,
                   unsigned int n,
                   unsigned int d,
                   std::vector<double> T_factors,
                    std::vector<double> F_factors,
                   std::vector<double> Ks )
{
  IMP::Particles children;
  IMP::Particles centroids;
  IMP::Particles images;
  IMP::Restraints R;
  Model* m=pf.get_model();
  int nlevels = ceil(log2(n));

  if (n==1)
    {
      // Create and return a singleton leaf:
      Particle* p =  pf.create("leaf %1%");
      return new TAMD_chain(p, centroids, images, R);
    }

  // Build <d> children recursively:
  {
    // factors and n particles in each child chain
    std::vector<double> T_factors1(T_factors.begin()+1, T_factors.end());
    std::vector<double> F_factors1(F_factors.begin()+1, F_factors.end());
    std::vector<double> Ks1(Ks.begin()+1, Ks.end());
    int per_child_base = n / d;
    int n_left = n;
    while(n_left > 0){
      int n_excess = n_left % d;
      int per_child = per_child_base + (n_excess > 0 ? 1 : 0);
      IMP::Pointer<TAMD_chain> child =
        create_tamd_chain(pf, per_child, d, T_factors, F_factors1, Ks1);
      n_left -= child.centroids.size();
      // Accumulate all results
      children.push_back( child.root );
      centroids.insert( centroids.end(),
                       child.centroids.begin(), child.centroids.end());
      images.insert ( images.end(),
                      child.images.begin(), child.images.end() );
      R.insert( R.end(),
                child.R.begin(), child.R.end());
    }
  }

  // Build centroid of <d> children
  Particle* pc;
  {
    std::ostringstream oss;  oss << "centroid " << nlevels;
    pc = new Particle( m, oss.str() );
    centroids.push_back( pc );
    core::Typed::setup_particle(ret_chain->root, pf.type_);
    core::XYZR pc_xyzr = core::XYZR::setup_particle(pc);
    pc_xyzr.set_radius(2); // TODO: something smarter?
    pc_xyzr.set_coordinates_are_optimized(true); // TODO: is needed or dangerous
                                                 // - for BD to evaluate it?
    atom::Diffusion pc_Diffusion = atom::Diffusion::setup_particle(pc);
    core::Hierarchy pc_h = core::Hierarchy::setup_particle(pc);
    for(unsigned int i = 0; i < children.size(); i++) {
      pc_h.add_child( IMP::core::Hierarchy(children[i]) );
    }
    base::Pointer<core::ChildrenRefiner> refiner =
      new core::ChildrenRefiner( core::Hierarchy::get_default_traits() );
    atom::CenterOfMass::setup_particle(pc, refiner.get());
    m->update(); // update now center of mass from children
  }

  // Build TAMD image of centroid + spring restraint
  {
    std::ostringstream oss;  oss << "TAMD " << nlevels << ".%1%";
    Particle* ptamd = create_tamd_image(  pc,
                                          oss.str(),
                                          T_factors[0],
                                          F_factors[0]);
    images.push_back(ptamd);
    base::Pointer<core::HarmonicDistancePairScore> spring=
      new core::HarmonicDistancePairScore(0, Ks[0]);
    R.push_back( new core::PairRestraint
                 ( spring, ParticlePair(pc, ptamd), oss.str() ) );
  }

  return new TAMD_chain(pc,centroids, images, R);
}



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
  IMP::core::Hierarchy::setup_particle( p_ret );
  IMP::atom::Diffusion::setup_particle( p_ret ); // diffusion coefficient?!
  IMP::atom::TAMDParticle::setup_particle( p_ret, T_factor, F_factor);
  IMP::atom::Mass::setup_particle( p_ret, atom::Mass(p_ref).get_mass() );
  core::Typed::setup_particle( p_ret, core::Typed(p_ref).get_type() );
  return p_ret;
}




IMPNPCTRANSPORT_END_INTERNAL_NAMESPACE
