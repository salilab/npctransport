/**
 *  ile creating_particles.cpp
 *  rief description.
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/creating_particles.h>
#include <IMP/core/XYZR.h>
#include <IMP/display/Colored.h>
#include <IMP/atom/Hierarchy.h>
#include <IMP/atom/Mass.h>
#include <IMP/core/rigid_bodies.h>
#include <IMP/atom/Diffusion.h>

#include "boost/tuple/tuple.hpp"


IMPNPCTRANSPORT_BEGIN_NAMESPACE
Particle* create_singleton_particle(SimulationData *sd, double radius,
                          double angular_D_factor,
                          double D_factor, display::Color c,
                          core::ParticleType type,
                          std::string name) {
  core::XYZR pc= core::XYZR::setup_particle(new Particle(sd->get_m()));
  pc->set_name(name);
  pc.set_radius(radius);
  pc.set_coordinates(algebra::get_zero_vector_d<3>());
  display::Colored::setup_particle(pc, c);
  core::Typed::setup_particle(pc, type);
  core::RigidBody rb
      =core::RigidBody::setup_particle(pc, algebra::ReferenceFrame3D());
  atom::RigidBodyDiffusion diff= atom::RigidBodyDiffusion::setup_particle(rb);
  double rdo=diff.get_rotational_diffusion_coefficient();
  diff.set_rotational_diffusion_coefficient(angular_D_factor*D_factor*rdo);
  diff.set_d(D_factor*diff.get_d());
  // rb.set_coordinates(IMP.algebra.get_random_vector_in(bb));
  rb.set_coordinates_are_optimized(true);
  pc->add_attribute(get_simulation_data_key(), sd);
  atom::Mass::setup_particle(pc, 1);
  return pc;
}

namespace {
/** Restraint the passed particles to be connected in a chain. The distance
    between consecutive particles is length_factor*the sum of the radii.
    // TODO: this documentation seems obsolete?

    Note, this assumes that all such chains will be disjoint and so you can
    use the container::ExclusiveConsecutivePairFilter if you want to filter
    out all pairs of particles connected by such chain restraints.

    The restraint is not added to the model.
*/
inline Restraint* create_chain_restraint(const ParticlesTemp &ps,
                                         LinearWellPairScore *pps,
                                         std::string name) {
  IMP_USAGE_CHECK(!ps.empty(), "No Particles passed.");

  // Exclusive means that the particles will be in no other
  // ConsecutivePairContainer
  // this assumption accelerates certain computations
  IMP_NEW(container::ExclusiveConsecutivePairContainer, cpc,
          (ps, name+" consecutive pairs"));
  Pointer<Restraint> r= container::create_restraint(pps, cpc.get(),
                                                    "chain restraint %1%");
  /\/ make sure it is not freed
  return r.release();
}
}

Particle*
create_chain(SimulationData *sd, int n, double radius,
             double angular_D_factor, double D_factor,
             LinearWellPairScore *ps,
             display::Color c,
             core::ParticleType t, std::string name) {
  ParticlesTemp ret;
  for ( int i=0; i< n; ++i) {
    ret.push_back(create_singleton_particle(sd, radius, angular_D_factor,
                                  D_factor, c, t, name));
  }
  Pointer<Restraint> cr=create_chain_restraint(ret,
                                               ps,
                                               name+"chain restraint");
  sd->add_chain_restraint(cr);
  atom::Hierarchy root
      = atom::Hierarchy::setup_particle(new Particle(sd->get_m()),
                                        ret);
  root->set_name(name);
  return root;
}



  /**        '''
        Create a TAMD image of centroid particle p

        p - reference particle to be tied by spring
        name - particle name
        T_factor - temeprature factor
        F_factor - friction factor

        returns the TAMD image particle
        '''
*/
Particle* create_tamd_image( Particle* p,
                             std::string name,
                             double T_factor,
                             double F_factor){
  //TAMD image of centroid
  Model* m = p->get_model();
  Particle* pstar = new IMP::Particle(m, name);
  IMP::core::XYZR pstarD = IMP::core::XYZR::setup_particle( pstar );
  IMP::core::XYZE pD = IMP::core::XYZR(p);
  pstarD.set_coordinates( pD.get_coordinates() );
  pstarD.set_radius( pD.get_radius() );
  pstarD.set_coordinates_are_optimized( True );
  IMP::core::Hierarchy pstarH = IMP::core::Hierarchy::setup_particle( pstar );
  IMP::atom::Diffusion.setup_particle( pstar ); // diffusion coefficient?!
  IMP::atom::TAMDParticle.setup_particle( pstar, T_factor, F_factor);
  IMP::atom::Mass pMass = IMP::atom::Mass( p );
  IMP::atom::Mass::setup_particle( pstar, pMass.get_mass() );
  return pstar;
}


/*
  Create a TAMD hierarchy of nlevels depth with a core for
  each d centroids in a lower level, with a real centroid and
  restrained centroid realization

  Params:
  -------
  m -         Model
  pf -        A factory for producing singleton particles
              (the leaves of the chain)
  nlevels -   Number of levels in the hierarchy. If 1 -
               return signleton particles
  d -         The out degree of each non-leaf node (# of children)
  T_factors - A list of length nlevels-1 with temeprature at each level
               from top to bottom
  G_factors - A list of length nlevels-1 with friction factor (G for gamma)
               at each level from top to bottom
  Ks -        Spring constants at each level between a particle and its TAMD
               image (a list of length nlevels-1)
        */
boost::tuple<
XXX _create_tamd_chain( Model* m,
                        ParticleFactory pf,
                        unsigned int nlevels,
                        unsinged int d,
                        std::vector<double> T_factors,
                        std::vector<double> F_factors,
                        std::vector<double> Ks )
{
  // build children hierarchies recursively first
  IMP::Particles children;
  IMP::Particles centroids;
  IMP::Particles images;
  IMP::Restraints R;  //  all restraints with image particles
  for (unsigned int i = 0; i < d; i++) {
    if (nlevels==1) {
      Particle* p =  pf.create("leaf %1%");
      Particles child_centroids;
      Particles *child_images;
      IMP::Restraints  child_R;
      return std::make_tuple<p, child_centroids, child_images, child_R>;
        }
    else {
                [p, child_centroids, child_images, child_R] \
                    = self._create_tamd_hierarchy(m,
                                                  nlevels - 1,
                                                  d,
                                                  T_factors[1:],
                                                  F_factors[1:],
                                                  Ks[1:])
                  }
            children.append(p)
            centroids = centroids + child_centroids
            images = images + child_images
            R = R + child_R
              }

        # Build centroid of lower hierarchies
        p = IMP.kernel.Particle(m, "centroid %d" % nlevels)
        centroids.append(p)
        pH = IMP.core.Hierarchy.setup_particle(p)
        pD = IMP.core.XYZR.setup_particle(p)
        pD.set_radius(2) # math.sqrt(math.pow(d, nlevels))) # TODO: something smarter?
#        pMass = IMP.atom.Mass.setup_particle(p, 1) # math.pow(d, nlevels)) # NOT RELEVANT FOR CENTER OF MASS
        pD.set_coordinates_are_optimized(True) # TODO: very unclear if this is needed or dangerous - to get BD to evaluate on it
        pDiffusion = IMP.atom.Diffusion.setup_particle(p)
        for child in children:
            pH.add_child( IMP.core.Hierarchy( child ) )
        print pH.get_children()
        refiner = IMP.core.ChildrenRefiner(
            IMP.core.Hierarchy.get_default_traits())
        IMP.atom.CenterOfMass.setup_particle(p, refiner)
        m.update() # so center of mass is up to date

        # Add a TAMD image of centroid with a spring
        pstar = self._create_tamd_image(  p,
                                          "TAMD image %d",
                                          T_factors[0],
                                          F_factors[0] )
        images.append(pstar)
        spring = IMP.core.HarmonicDistancePairScore(0, Ks[0])
        r=IMP.core.PairRestraint(spring, (p, pstar))
        R.append(r)

        return p, centroids, images, R

    def _get_ordered_leaves(self, h):
        '''
        get leave particles by certain order that is guaranteed to cluster
        leave particles with a common ancestor together, for hierarchy h
        '''
        if(h.get_number_of_children() == 0):
            return [h.get_particle()]
        leaves = []
        for child in h.get_children():
            leaves = leaves + self._get_ordered_leaves(child)
        return leaves


IMPNPCTRANSPORT_END_NAMESPACE
