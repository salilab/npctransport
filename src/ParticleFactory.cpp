/**
 * \file ParticleFactory.cpp
 * \brief create singleton particles for npc simulations
 *
 *  Copyright 2007-2018 IMP Inventors. All rights reserved.
 */

#include <IMP/npctransport/ParticleFactory.h>
#include <IMP/core/Typed.h>
#include <IMP/core/XYZR.h>
#include <IMP/display/Colored.h>
#include <IMP/atom/Hierarchy.h>
#include <IMP/atom/Mass.h>
#include <IMP/core/rigid_bodies.h>
#include <IMP/atom/Diffusion.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE


/**
   create a particle associated with the model of sd, with params
   as specified in factory class
 */
Particle*
ParticleFactory::create(std::string name) {
  set_was_used(true);
  if(name==""){
    name = type_.get_string();
  }
  Particle* p = new Particle(get_model(), name);

  core::XYZR p_xyzr= core::XYZR::setup_particle(p);
  p_xyzr.set_radius(radius_);
  p_xyzr.set_coordinates(algebra::get_zero_vector_d<3>());

  display::Colored::setup_particle(p, color_);

  core::Typed::setup_particle(p, type_);

  core::RigidBody p_rb
      =core::RigidBody::setup_particle(p, algebra::ReferenceFrame3D());
  p_rb.set_coordinates_are_optimized(true);

  if(D_factor_ * angular_D_factor_ > 0.0) {
    atom::RigidBodyDiffusion p_rbd= atom::RigidBodyDiffusion::setup_particle(p);
    p_rbd.set_diffusion_coefficient
      ( D_factor_ * p_rbd.get_diffusion_coefficient() );
    double rdc=p_rbd.get_rotational_diffusion_coefficient();
    p_rbd.set_rotational_diffusion_coefficient
      ( angular_D_factor_ * D_factor_ * rdc );
  } else if ( D_factor_ > 0.0 ){
    atom::Diffusion p_d= atom::Diffusion::setup_particle(p);
    p_d.set_diffusion_coefficient
      ( D_factor_ * p_d.get_diffusion_coefficient() );
  }

  atom::Mass::setup_particle(p, 1);

  p->add_attribute(get_simulation_data_key(), sd_);

  return p;
}

IMPNPCTRANSPORT_END_NAMESPACE
