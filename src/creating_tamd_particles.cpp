/**
 *  \file creating_particles.cpp
 *  \brief description.
 *
 *  Copyright 2007-2014 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/creating_tamd_particles.h>
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
#include <cstdio>
#include <sstream>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

Particle *create_tamd_fg_chain
( SimulationData *sd,
  const ::npctransport_proto::Assignment_FGAssignment &fg_data,
  display::Color c )
{
  // create particles

  core::ParticleType type(fg_data.type());
  int n = fg_data.number_of_beads().value();
  double radius = fg_data.radius().value();
  double D_factor = fg_data.d_factor().value();
  double angular_D_factor = sd->get_angular_d_factor();
  ParticleFactory pf(sd, radius, D_factor, angular_D_factor, c, type);
  ParticlesTemp particles;
  for (int i = 0; i < n; ++i) {
    particles.push_back( pf.create() );
  }
  // put in hierarchy
  atom::Hierarchy root =
    atom::Hierarchy::setup_particle( new Particle(sd->get_model() ),
                                     particles );
  std::string name = type.get_string();
  root->set_name(name);
  // add chain backbone restraint
  double rlf = fg_data.rest_length_factor().value();
  sd->get_scoring()->add_chain_restraint
    ( root, rlf, name + "chain restraint" );

  return root;
}

namespace {

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
  Model* m = p_ref->get_model();
  Particle* p_ret = new IMP::Particle(m, name);
  IMP::core::XYZR p_ret_xyzr = IMP::core::XYZR::setup_particle( p_ret );
  IMP::core::XYZR p_ref_xyzr = IMP::core::XYZR(p_ref);
  p_ret_xyzr.set_coordinates( p_ref_xyzr.get_coordinates() );
  p_ret_xyzr.set_radius( p_ref_xyzr.get_radius() );
  p_ret_xyzr.set_coordinates_are_optimized( true );
  IMP::core::Hierarchy p_ret_h = IMP::core::Hierarchy::setup_particle( p_ret );
  IMP::atom::Diffusion::setup_particle( p_ret ); // diffusion coefficient?!
  IMP::atom::TAMDParticle::setup_particle( p_ret, T_factor, F_factor);
  IMP::atom::Mass p_ref_mass = IMP::atom::Mass( p_ref );
  IMP::atom::Mass::setup_particle( p_ret, p_ref_mass.get_mass() );
  return p_ret;
}
}

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
                   std::vector<double> Ks )
{
  IMP::Particles children;
  IMP::Particles centroids;
  IMP::Particles images;
  IMP::Restraints R;

  if (nlevels==0)
    {
      // Create and return a singleton leaf:
      Particle* p =  pf.create("leaf %1%");
      return TAMD_chain(p, centroids, images, R);
                             // IMP::Particles(),
                             // IMP::Particles(),
                             // IMP::Restraints()>;
    }

  // Build <d> children recursively:
  {
    std::vector<double> T_factors1(T_factors.begin()+1, T_factors.end());
    std::vector<double> F_factors1(F_factors.begin()+1, F_factors.end());
    std::vector<double> Ks1(Ks.begin()+1, Ks.end());
    for (unsigned int i = 0; i < d; i++) {
      TAMD_chain tc = create_tamd_chain(m, pf, nlevels - 1, d,
                                        T_factors, F_factors1, Ks1);
      // Accumulate all results
      children.push_back( tc.p );
      std::copy ( tc.centroids.begin(),
                  tc.centroids.end(),
                  back_inserter(centroids) );
      std::copy ( tc.images.begin(),
                  tc.images.end(),
                  back_inserter(images) );
      std::copy ( tc.R.begin(),
                  tc.R.end(),
                  back_inserter(R) );
    }
  }

  // Build centroid of <d> children
  Particle* pc;
  {
    std::ostringstream oss;  oss << "centroid " << nlevels;
    pc = new Particle( m, oss.str() );
    centroids.push_back( pc );
    core::XYZR pc_xyzr = core::XYZR::setup_particle(pc);
    pc_xyzr.set_radius(2); // TODO: something smarter?
    pc_xyzr.set_coordinates_are_optimized(true); // TODO: is needed or dangerous
                                                 // - for BD to evaluate it?
    // IMP::atom::Mass.setup_particle(pc, 1) // NOT RELEVANT FOR CENTER OF MASS
    atom::Diffusion pc_Diffusion = atom::Diffusion::setup_particle(pc);
    core::Hierarchy pc_h = core::Hierarchy::setup_particle(pc);
    for(unsigned int i = 0; i < children.size(); i++) {
      pc_h.add_child( IMP::core::Hierarchy(children[i]) );
    }
    //    std::cout << pc_h.get_children() << std::endl;
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

  return TAMD_chain(pc,centroids, images, R);
}

    // def _get_ordered_leaves(self, h):
    //     '''
    //     get leave particles by certain order that is guaranteed to cluster
    //     leave particles with a common ancestor together, for hierarchy h
    //     '''
    //     if(h.get_number_of_children() == 0):
    //         return [h.get_particle()]
    //     leaves = []
    //     for child in h.get_children():
    //         leaves = leaves + self._get_ordered_leaves(child)
    //     return leaves


IMPNPCTRANSPORT_END_NAMESPACE
