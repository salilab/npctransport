/**
 *  \file creating_particles.cpp
 *  \brief description.
 *
 *  Copyright 2007-2015 IMP Inventors. All rights reserved.
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
   ( particles, rlf, name + "chain restraint" );

  return root;
}

IMPNPCTRANSPORT_END_NAMESPACE
