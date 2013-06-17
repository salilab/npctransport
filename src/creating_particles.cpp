/**
 *  ile creating_particles.cpp
 *  rief description.
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/creating_particles.h>
#include <IMP/container/ConsecutivePairContainer.h>
#include <IMP/core/XYZR.h>
#include <IMP/display/Colored.h>
#include <IMP/atom/Hierarchy.h>
#include <IMP/atom/Mass.h>
#include <IMP/core/rigid_bodies.h>
#include <IMP/atom/Diffusion.h>
//#include <IMP/example/creating_restraints.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE
Particle *create_particle(SimulationData *sd, double radius,
                          double angular_D_factor, double D_factor,
                          display::Color c, core::ParticleType type,
                          std::string name) {
  core::XYZR pc = core::XYZR::setup_particle(new Particle(sd->get_m()));
  pc->set_name(name);
  pc.set_radius(radius);
  pc.set_coordinates(algebra::get_zero_vector_d<3>());
  display::Colored::setup_particle(pc, c);
  core::Typed::setup_particle(pc, type);
  core::RigidBody rb =
      core::RigidBody::setup_particle(pc, algebra::ReferenceFrame3D());
  atom::RigidBodyDiffusion diff = atom::RigidBodyDiffusion::setup_particle(rb);
  double rdo = diff.get_rotational_diffusion_coefficient();
  diff.set_rotational_diffusion_coefficient(angular_D_factor * D_factor * rdo);
  diff.set_d(D_factor * diff.get_d());
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
inline Restraint *create_chain_restraint(const ParticlesTemp &ps,
                                         LinearWellPairScore *pps,
                                         std::string name) {
  IMP_USAGE_CHECK(!ps.empty(), "No Particles passed.");

  // Exclusive means that the particles will be in no other
  // ConsecutivePairContainer
  // this assumption accelerates certain computations
  IMP_NEW(container::ExclusiveConsecutivePairContainer, cpc,
          (ps, name + " consecutive pairs"));
  base::Pointer<Restraint> r =
      container::create_restraint(pps, cpc.get(), "chain restraint %1%");
  // make sure it is not freed
  return r.release();
}
}

Particle *create_chain(SimulationData *sd, int n, double radius,
                       double angular_D_factor, double D_factor,
                       LinearWellPairScore *ps, display::Color c,
                       core::ParticleType t, std::string name) {
  ParticlesTemp ret;
  for (int i = 0; i < n; ++i) {
    ret.push_back(
        create_particle(sd, radius, angular_D_factor, D_factor, c, t, name));
  }
  base::Pointer<Restraint> cr =
      create_chain_restraint(ret, ps, name + "chain restraint");
  sd->add_chain_restraint(cr);
  atom::Hierarchy root =
      atom::Hierarchy::setup_particle(new Particle(sd->get_m()), ret);
  root->set_name(name);
  return root;
}
IMPNPCTRANSPORT_END_NAMESPACE
