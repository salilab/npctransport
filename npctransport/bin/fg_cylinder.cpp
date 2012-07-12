/**
 * \file fg_array.cpp
 * \brief Simulate an fg and a kap interacting
 * Copyright 2007-2012 IMP Inventors. All rights reserved.
 */

#define IMP_NPC_MAIN
#include <IMP/npctransport/main.h>

#include <RMF/utility.h>


int main(int argc, char *argv[]) {
  using namespace IMP;
  using namespace IMP::npctransport;
  IMP_NPC_STARTUP(sd);
  RMF:: set_show_hdf5_errors(true);
  //  sd->add_interaction(type_of_float[0], type_of_fg[0]);
  //  sd->add_interaction(type_of_float[1], type_of_fg[0]);

  //place fg anchors in grid and don't optimize them
  algebra::Cylinder3D cyl= sd->get_cylinder();

  atom::Hierarchy root= sd->get_root();
  atom::Hierarchies chains =get_fg_chains(root);
  // create a set of random sites (for now)
  algebra::Vector3Ds sites=get_grid_surface_cover(cyl, sqrt(chains.size())+1,
                                                  sqrt(chains.size())+1);
  std::cout << IMP::base::Showable(sites) << std::endl;
  // anchor each fg chain to the (x,y) site (only by x,y coords,
  // all anchored to the same z plane)
  for (unsigned int i=0; i< chains.size(); ++i) {
    atom::Hierarchy r(chains[i]);
    core::XYZ d(r.get_child(0));
    d.set_coordinates(sites[i]);
    d.set_coordinates_are_optimized(false);
  }
  IMP_NPC_LOOP(sd, ParticlePairsTemp());
  return 0;
}
