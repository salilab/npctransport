/**
 *  \file simulation_data.cpp
 *  \brief description.
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/npctransport/particle_types.h>

IMPNPCTRANSPORT_BEGIN_NAMESPACE

// TODO: the specific type names ("fg0", etc.)  should all be user parameters
core::ParticleType type_of_fga[] = {
  core::ParticleType("fg0"), core::ParticleType("fg1"),
  core::ParticleType("fg2"), core::ParticleType("fg3"),
  core::ParticleType("fg4"), core::ParticleType("fg5"),
  core::ParticleType("fg6"), core::ParticleType("fg7"),
  core::ParticleType("fg8"), core::ParticleType("fg9"),
  core::ParticleType("fg10"), core::ParticleType("fg11"),
  core::ParticleType("fg12"), core::ParticleType("fg13"),
  core::ParticleType("fg14"), core::ParticleType("fg15")
};

core::ParticleType type_of_floata[] = {core::ParticleType("kap"),
                                       core::ParticleType("crap0"),
                                       core::ParticleType("crap1")};
core::ParticleTypes type_of_fg(type_of_fga,
                               type_of_fga + sizeof(type_of_fga) /
                                                 sizeof(core::ParticleType));
core::ParticleTypes type_of_float(type_of_floata,
                                  type_of_floata +
                                      sizeof(type_of_floata) /
                                          sizeof(core::ParticleType));

atom::Hierarchies get_fg_chains(atom::Hierarchy root) {
  atom::Hierarchies ret;
  // I. return root itself if the type of its first direct child
  // is contained in [type_of_fg]
  if (root.get_number_of_children() > 0) {
    atom::Hierarchy c = root.get_child(0);
    if (core::Typed::get_is_setup(c)) {
      core::ParticleType t = core::Typed(c).get_type();
      if (std::find(type_of_fg.begin(), type_of_fg.end(), t) !=
          type_of_fg.end()) {
        return atom::Hierarchies(1, root);
      }
    }
  }
  // II. Otherwise, recurse on all of root's children
  for (unsigned int i = 0; i < root.get_number_of_children(); ++i) {
    ret += get_fg_chains(root.get_child(i));
  }
  return ret;
}

IMPNPCTRANSPORT_END_NAMESPACE
