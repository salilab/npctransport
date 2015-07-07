/**
# * \file fg_simulation.cpp
# * \brief Simulate an fg and a kap interacting
#
# * Copyright 2007-2012 IMP Inventors. All rights reserved.
# */

//#include <IMP/npctransport/npctransport_config.h>
#include <IMP/npctransport/main.h>
#include <IMP/npctransport/FGChain.h>
#include <IMP/atom/Hierarchy.h>
#include <IMP/algebra/Vector2D.h>
#include <IMP/algebra/Vector3D.h>
#include <IMP/algebra/BoundingBoxD.h>
#include <IMP/display/Colored.h>
#include <IMP/display/Color.h>
#include <IMP/container/SingletonsRestraint.h>
#include <IMP/CreateLogContext.h>
#include <IMP/exception.h>
#include <IMP/enums.h>
#include <IMP/npctransport.h>
#include <IMP/base_types.h>
#include <IMP/Pointer.h>
#include <IMP/base_statistics.h>
#include <IMP/Restraint.h>
#include <IMP/SingletonScore.h>
#include <IMP/core/Typed.h>
#include <numeric>
#include <cmath>
#include <iostream>
#include <IMP/npctransport/internal/boost_main.h>

// use the example code for now to work bugs out of it
//#include <IMP/example/creating_restraints.h>
//#include <IMP/example/randomizing.h>
//#include <IMP/example/counting.h>
//#include <IMP/example/optimizing.h>

namespace {
/**
    anchors all the fgs to a planar surface
    at the edge of the simulation bounding box

    @param sd Simulation data
*/
void set_fg_grid(IMP::npctransport::SimulationData& sd) {
  using namespace IMP;
  using namespace IMP::algebra;

  //  get bottom surface of the simulation data bounding box
  Vector2D lower_corner_XY(sd.get_box().get_corner(0)[0],
                           sd.get_box().get_corner(0)[1]);
  Vector2D upper_corner_XY(sd.get_box().get_corner(1)[0],
                           sd.get_box().get_corner(1)[1]);
  algebra::BoundingBox2D surface(lower_corner_XY, upper_corner_XY);
  // get fg
  atom::Hierarchies chain_roots = sd.get_fg_chain_roots();
  // anchor fgs to surface,
  // for now using random non-overlapping points
  // create a set of random sites (for now)
  Pointer<npctransport::FGChain> first_chain =
    npctransport::get_fg_chain(chain_roots[0]);
  double r = core::XYZR(first_chain->get_bead(0)).get_radius();
  Vector2Ds sites;
  std::cout << IMP::Showable(sites) << std::endl;
  while (sites.size() < chain_roots.size()) {
    Vector2D cur = get_random_vector_in(surface);
    bool bad = false;
    for (unsigned int i = 0; i < sites.size(); ++i) {
      // 2*r non-overlapping
      if (get_distance(sites[i], cur) < 2 * r) {
        bad = true;
        break;
      }
    }
    if (!bad) {
      sites.push_back(cur);
      std::cout << "Site # " << sites.size() << " is " << cur << std::endl;
    }
  }
  // anchor each fg chain to a site
  for (unsigned int i = 0; i < chain_roots.size(); ++i) {
    Pointer<npctransport::FGChain> chain =
      npctransport::get_fg_chain((chain_roots[i]));
    core::XYZ d(chain->get_bead(0));
    d.set_coordinates(
        Vector3D(sites[i][0], sites[i][1], sd.get_box().get_corner(0)[2]));
    d.set_coordinates_are_optimized(false);
    std::cout << "d = " << d << std::endl;
  }
}

/*
  color the different fgs in different colors

  @param chains  the SimulationData object
*/
void color_fgs(IMP::npctransport::SimulationData& sd) {
  using namespace IMP;
  using namespace IMP::npctransport;
  using IMP::display::Colored;

  atom::Hierarchies chain_roots(sd.get_fg_chain_roots());
  unsigned int n_chains = chain_roots.size();
  for (unsigned int i = 0; i < n_chains; i++) {
    display::Color color;
    // choose color
    if (n_chains <= 11) {
      color = display::get_display_color(i);
    } else {
      double f = i / (float)(n_chains - 1);  // spread in [0..1]
      color = display::get_jet_color(f);
    }
    // apply color
    atom::Hierarchies children = core::get_all_descendants(chain_roots[i]);
    for (unsigned int j = 0; j < children.size(); j++) {
      if (Colored::get_is_setup(children[j])) {
        Colored(children[j]).set_color(color);
      } else {
        Colored::setup_particle(children[j], color);
      }
    }
  }
}

/**
   anchors the FGs to the surface of the simulation bounding cylinder
   (= slab constraint)

   sd - the SimulationData object
   n_layers - number of fg nup layers
*/
void set_fgs_in_cylinder(IMP::npctransport::SimulationData& sd, int n_layers) {
  using namespace IMP;
  using atom::Hierarchy;
  using atom::Hierarchies;

  IMP::algebra::Cylinder3D cyl = sd.get_cylinder();
  Hierarchies chain_roots = sd.get_fg_chain_roots();
  // compute the relative radius in which particles would be positioned
  // TODO: we assume here that particle radius is smaller
  //       than the cylinder radius - verify in runtime?
  Particles chain_beads =
    npctransport::get_fg_chain(chain_roots[0])->get_beads();
  double particle_radius = IMP::core::XYZR(chain_beads[0]).get_radius();
  // compute fraction of particle from full cylinder radius
  double relative_r = (cyl.get_radius() - particle_radius) / cyl.get_radius();
  // compute vertical poisition along central axis, and inter-layer distance
  double h_bottom_layer;
  if (n_layers == 1)
    h_bottom_layer = 0.5;
  else
    h_bottom_layer = 0.0;
  double dLayers = 0.0;
  if (n_layers > 1) {
    dLayers = 1.0 / (n_layers - 1);
  }
  // calculate angle increments between adjacent fg nups in each layers
  int chains_per_layer = (int)(std::ceil(chain_roots.size() / (n_layers + 0.0)));
  double angle_increments = 2 * IMP::PI / chains_per_layer;
  // pin chains to each layer
  for (int layer_num = 0; layer_num < n_layers; layer_num++) {
    double relative_h = h_bottom_layer + layer_num * dLayers;
    for (int k = 0; k < chains_per_layer; k++) {
      unsigned int chain_num = layer_num * chains_per_layer + k;
      if (chain_num >= chain_roots.size())
        break;  // may happen if len(chains) does not divide by n_layers
      double angle = k * angle_increments;
      algebra::Vector3D new_anchor =
          cyl.get_inner_point_at(relative_h, relative_r, angle);
      Pointer<npctransport::FGChain> cur_chain =
        npctransport::get_fg_chain(chain_roots[chain_num]);
      core::XYZ d(cur_chain->get_bead(0));
      d.set_coordinates(new_anchor);
      d.set_coordinates_are_optimized(false);
      std::cout << "d = " << d << std::endl;
    }
  }
}

/**    returns all kap / crap particles in SimulationData */
IMP::ParticlesTemp get_kaps_and_craps(IMP::npctransport::SimulationData& sd) {
  using namespace IMP::npctransport;

  IMP::ParticlesTemp ret;
  ParticleTypeSet const &fts = sd.get_floater_types();
  for (ParticleTypeSet::const_iterator it = fts.begin();
       it != fts.end();
       it++)
    {
      std::cout << *it << std::endl;
      ret += sd.get_root_of_type(*it).get_children();
      std::cout << ret;
    }
  return ret;
}

IMP::Pointer<IMP::Restraint> get_exclude_from_channel_restraint(
    IMP::npctransport::SimulationData& sd) {
  using namespace IMP;

  bool both_sides = ! sd.get_are_floaters_on_one_slab_side();
  double top = (sd.get_slab_thickness() / 2) * 1.5;  // *1.5 to get some slack
  double LARGE = 10000000;
  double bottom = both_sides ? -top : (-top * LARGE);
  double k = 40.0;
  IMP_NEW(IMP::npctransport::ExcludeZRangeSingletonScore, score,
          (bottom, top, k));
  IMP::ParticlesTemp particles = get_kaps_and_craps(sd);
  IMP_NEW(container::SingletonsRestraint, sr,
          (score, particles, "ExcludeZRangeRestraint"));
  return sr;
}

  // print the first atoms of all the fgs in sd
  void print_fgs(IMP::npctransport::SimulationData& sd, IMP::LogLevel ll)
  {
    using namespace IMP;
    using atom::Hierarchy;
    using atom::Hierarchies;

    Hierarchies chains = sd.get_fg_chain_roots();
    for (unsigned int k = 0; k < chains.size(); k++) {
      Pointer<npctransport::FGChain> cur_chain =
        npctransport::get_fg_chain(chains[k]);
      core::XYZ d(cur_chain->get_bead(0));
      IMP_LOG(ll, "d # " << k << " = " << d << std::endl);
      IMP_LOG(ll, "is optimizable = "
              << d.get_coordinates_are_optimized()
              << std::endl);
    }
  }

boost::int64_t cylinder_nlayers = 0;
IMP::AddIntFlag cylinder_adder(
    "cylinder_nlayers", "anchor FG nups to a cylindrical pore with a slab of"
                        " dimensions that are specified in"
                        " the config file, with specified number of FG layers"
                        " (no cylinder anchoring if equals to 0)",
    &cylinder_nlayers);

bool surface_anchoring = false;
IMP::AddBoolFlag surface_adder("surface_anchoring",
                                     "anchor FG nups to bottom of cube",
                                     &surface_anchoring);
}


int main(int argc, char* argv[]) {
  using namespace IMP;

  // logging stuff:
  IMP::CreateLogContext main("main");
  // preparation::
  try {
    IMP::Pointer<npctransport::SimulationData> sd =
        npctransport::startup(argc, argv);
    //print_fgs(*sd, IMP::TERSE);
    if (surface_anchoring) {
      IMP_ALWAYS_CHECK(cylinder_nlayers == 0,
                       "surface anchoring and cylinder"
                       " anchoring flags are mutually exclusive",
                       IMP::ValueException);
      set_fg_grid(*sd);
    }
    if (cylinder_nlayers > 0) {
      set_fgs_in_cylinder(*sd, cylinder_nlayers);
      std::cout << "Numbner of cylinder layers = " << cylinder_nlayers
                << std::endl;
    }
    color_fgs(*sd);
    Restraints initialization_restraints;
    if (sd->get_scoring()->get_has_slab() &&
        sd->get_is_exclude_floaters_from_slab_initially()) {
      // if has slab, exclude from channel initially
      IMP::Pointer<IMP::Restraint> r =
          get_exclude_from_channel_restraint(*sd);
      initialization_restraints.push_back(r);
    }
    //IMP_LOG(PROGRESS, initialization_restraints << std::endl);
    npctransport::do_main_loop(sd, initialization_restraints);
    //print_fgs(*sd, TERSE);
    IMP::show_timings(std::cout);
  }
  catch (const IMP::Exception & e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return -1;
  }
  return 0;
}
