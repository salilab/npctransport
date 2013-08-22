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
#include "npctransport.pb.h"

//#include <IMP/example/creating_restraints.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE
Particle *create_particle(SimulationData *sd, double radius,
                          double angular_D_factor, double D_factor,
                          display::Color c, core::ParticleType type) {
  core::XYZR pc = core::XYZR::setup_particle(new Particle(sd->get_model()));
  std::string name = type.get_string();
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
  diff.set_diffusion_coefficient(D_factor * diff.get_diffusion_coefficient());
  // rb.set_coordinates(IMP.algebra.get_random_vector_in(bb));
  rb.set_coordinates_are_optimized(true);
  pc->add_attribute(get_simulation_data_key(), sd);
  atom::Mass::setup_particle(pc, 1);
  return pc;
}

Particle *create_fg_chain
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
 ParticlesTemp particles;
 for (int i = 0; i < n; ++i) {
   particles.push_back(
        create_particle(sd, radius, angular_D_factor, D_factor, c, type));
 }
 // put in hierarchy
  atom::Hierarchy root =
    atom::Hierarchy::setup_particle( new Particle(sd->get_model() ),
                                     particles );
  std::string name = type.get_string();
  root->set_name(name);
  // add chain backbone restraint
  double rlf = fg_data.rest_length_factor().value();
  double rest_length = 2 * radius * rlf;
  sd->get_scoring()->add_chain_restraint
    ( root, rest_length, name + "chain restraint" );

  return root;
}

IMPNPCTRANSPORT_END_NAMESPACE
