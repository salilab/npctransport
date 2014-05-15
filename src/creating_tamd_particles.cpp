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
#include <cmath>
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
  int d = 2; // outdegree in TAMD hierarchy, TODO: parametrize
  int n_levels = ceil(log(n)/log(d)); // log_d{n}
  std::vector<double> T_factors(n_levels); // temperature scaling
  std::vector<double> F_factors(n_levels); // friction scaling
  std::vector<double> Ks(n_levels); // TAMD spring constant
  for(int i=0; i < n_levels; i++) { // i ~ increasing depth from root
    int level = n_levels - i; // level above leaves
    T_factors[i] = 3 * pow(2,level-1);
    F_factors[i] = 15 * pow(3,level-1);
    Ks[i] = 10;
  }
  TAMD_chain tc = create_tamd_chain(pf, n, d,
                                    T_factors, F_factors, Ks);

  // add chain backbone restraint
  double rlf = fg_data.rest_length_factor().value();
  Particles L = get_ordered_tamd_leaves( core::Hierarchy(tc.root) );
  sd->get_scoring()->add_chain_restraint
    ( L, rlf, type.get_string() + "chain restraint" );

  return tc.root;
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
  @param n       Number of particles in chain (=leaves). If 1, return singleton.
  @param d       maximal out degree of each non-leaf node (# of children)
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
      return TAMD_chain(p, centroids, images, R);
                             // IMP::Particles(),
                             // IMP::Particles(),
                             // IMP::Restraints()>;
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
      TAMD_chain child = create_tamd_chain(pf, per_child, d,
                                        T_factors, F_factors1, Ks1);
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

/**
   get leave particles by certain order that is guaranteed to cluster
   leave particles with a common ancestor together, for hierarchy h
*/
Particles get_ordered_tamd_leaves(core::Hierarchy root) {
  int n = root.get_number_of_children();
  if(n == 0){
    return Particles(1, root.get_particle());
  }
  Particles leaves;
  for(int i = 0 ; i < n; i++) {
    Particles sub_leaves = get_ordered_tamd_leaves(root.get_child(i));
    leaves.insert(leaves.end(),
                  sub_leaves.begin(), sub_leaves.end());
  }
  return leaves;
}


IMPNPCTRANSPORT_END_NAMESPACE
