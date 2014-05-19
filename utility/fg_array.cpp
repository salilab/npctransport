/**
 * \file fg_array.cpp
 * \brief Simulate an fg and a kap interacting
 * Copyright 2007-2012 IMP Inventors. All rights reserved.
 */

#define IMP_NPC_MAIN
#include <IMP/npctransport/main.h>
#include <IMP/ParticleTuple.h>
#include <IMP/npctransport/FGChain.h>
#include <IMP/base/Pointer.h>
#include <RMF/utility.h>

namespace {
int do_it(IMP::base::Pointer<IMP::npctransport::SimulationData> sd) {
  using namespace IMP;
  using namespace IMP::npctransport;
  using namespace IMP::algebra;
  //  sd->add_interaction(type_of_float[0], type_of_fg[0]);
  //  sd->add_interaction(type_of_float[1], type_of_fg[0]);

  // place fg anchors in grid and don't optimize them
  IMP::algebra::BoundingBox2D base(
      IMP::algebra::Vector2D(sd->get_box().get_corner(0)[0],
                             sd->get_box().get_corner(0)[1]),
      IMP::algebra::Vector2D(sd->get_box().get_corner(1)[0],
                             sd->get_box().get_corner(1)[1]));
  IMP::atom::Hierarchy root = sd->get_root();
  IMP::atom::Hierarchies chain_roots = sd->get_fg_chain_roots();
  // create a set of random sites (for now)
  IMP::algebra::Vector2Ds sites;
  std::cout << IMP::base::Showable(sites) << std::endl;
  base::Pointer<FGChain> chain = get_fg_chain(chain_roots[0]);
  double r = IMP::core::XYZR(chain->beads[0]).get_radius();
  std::cout << "Base is " << base << std::endl;
  do {
    // add a site that is not too close to an existing sites (at least 2*r)
    // - a site for every cain
    IMP::algebra::Vector2D cur = (get_random_vector_in(base));
    bool bad = false;
    for (unsigned int i = 0; i < sites.size(); ++i) {
      if (get_distance(sites[i], cur) < 2 * r) {
        bad = true;
        break;
      }
    }
    if (!bad) {
      sites.push_back(cur);
    }
  } while (sites.size() < chain_roots.size());
  // anchor each fg chain to the (x,y) site (only by x,y coords,
  // all anchored to the same z plane)
  for (unsigned int i = 0; i < chain_roots.size(); ++i) {
    base::Pointer<FGChain> chain = get_fg_chain(chain_roots[i]);
    IMP::core::XYZ d(chain->beads[0]);
    d.set_coordinates(
        Vector3D(sites[i][0], sites[i][1], sd->get_box().get_corner(0)[2]));
    d.set_coordinates_are_optimized(false);
  }
  IMP::npctransport::do_main_loop(sd, IMP::RestraintsTemp());
  return 0;
}
}

int main(int argc, char *argv[]) {
  int ret;
  IMP::base::Pointer<IMP::npctransport::SimulationData> sd =
      IMP::npctransport::startup(argc, argv);

  ret = do_it(sd);

  return ret;
}
