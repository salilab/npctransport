/**
# * \file fg_simulation.cpp
# * \brief Simulate an fg and a kap interacting
#
# * Copyright 2007-2012 IMP Inventors. All rights reserved.
# */

#include <IMP/npctransport/npctransport_config.h>
#include <IMP/algebra.h>
#include <IMP/core.h>
#include <IMP/atom.h>
#include <IMP/display.h>
#include <IMP/rmf.h>
#include <IMP.h>
#include <RMF/utility.h>
#include <IMP/container.h>
#include <IMP/base/CreateLogContext.h>
#include <IMP/npctransport.h>
#include <IMP/benchmark/Profiler.h>
#include <numeric>
#include <cmath>
#include <iostream>

// use the example code for now to work bugs out of it
#include <IMP/example/creating_restraints.h>
#include <IMP/example/randomizing.h>
#include <IMP/example/counting.h>
#include <IMP/example/optimizing.h>

#define IMP_NPC_MAIN
#ifdef IMP_NPC_GOOGLE
#include <IMP/npctransport/internal/google_main.h>
#else
#include <IMP/npctransport/internal/boost_main.h>
#endif


IMP_NPC_PARAMETER_INT(work_unit, -1, "The work unit");
IMP_NPC_PARAMETER_STRING(configuration, "configuration.pb",
                         "input configuration file in protobuf format"
                         " [default: %default]");
IMP_NPC_PARAMETER_STRING(assignments, "assignments.pb",
                         "output assignments file in protobuf format,"
                         " recording the assignment being executed"
                         " [default: %default]");
IMP_NPC_PARAMETER_STRING(statistics, "statistics.pb",
                         "output statistics file in protobuf format"
                         " [default: %default]");
IMP_NPC_PARAMETER_STRING(final_configuration, "final.pym",
                         "output final configuration file"
                         " [default: %default]");
IMP_NPC_PARAMETER_STRING(rmf_file, "output.rmf",
                         "RMF file for recording the simulation progress"
                         " [default: %default]");
IMP_NPC_PARAMETER_BOOL(profile, false,
                       "Whether to turn on profiling for the first run");
IMP_NPC_PARAMETER_BOOL(quick, false,
                       "Reduce all steps to the minimum");
IMP_NPC_PARAMETER_BOOL(show_steps, false,
                       "Show the steps for each modified variable");
IMP_NPC_PARAMETER_BOOL(show_number_of_work_units, false,
                       "Show the number of work units");
IMP_NPC_PARAMETER_BOOL(cylinder_anchoring, false,
                       "anchor FG nups to a cylinder specified in"
                       "the config file");


/**
    anchors all the fgs to a planar surface
    at the edge of the simulation bounding box

    @param sd Simulation data
*/
void set_fg_grid(IMP::npctransport::SimulationData& sd )
{
  using namespace IMP;
  using namespace IMP::algebra;

  //  get bottom surface of the simulation data bounding box
  Vector2D lower_corner_XY
    (sd.get_box().get_corner(0)[0],
     sd.get_box().get_corner(0)[1]) ;
  Vector2D upper_corner_XY
    (sd.get_box().get_corner(1)[0],
     sd.get_box().get_corner(1)[1]) ;
  BoundingBox2D surface
    ( lower_corner_XY, upper_corner_XY ) ;
  // get fg
  atom::Hierarchy root= sd.get_root();
  atom::Hierarchies chains = IMP::npctransport::get_fg_chains(root);
  // anchor fgs to surface,
  // for now using random non-overlapping points
  // create a set of random sites (for now)
  double r= core::XYZR(chains[0].get_child(0)).get_radius();
  Vector2Ds sites;
  std::cout << base::Showable(sites) << std::endl;
  while ( sites.size() < chains.size() )
    {
      Vector2D cur = get_random_vector_in(surface);
      bool bad = false;
      for (unsigned int i=0; i< sites.size(); ++i) {
        // 2*r non-overlapping
        if (get_distance(sites[i], cur) < 2*r) {
          bad=true;
          break;
        }
      }
      if (!bad) {
        sites.push_back(cur);
        std::cout << "Site # " << sites.size()
                  << " is " << cur << std::endl;
      }
    }
  // anchor each fg chain to a site
  for (unsigned int i=0; i< chains.size(); ++i) {
    atom::Hierarchy r(chains[i]);
    core::XYZ d(r.get_child(0));
    d.set_coordinates(Vector3D(sites[i][0], sites[i][1],
                               sd.get_box().get_corner(0)[2]));
    d.set_coordinates_are_optimized( false );
    std::cout << "d = " << d << std::endl;
  }
}


/*
  color the different fgs in different colors

  @param chains  the SimulationData object
*/
void color_fgs( IMP::npctransport::SimulationData& sd ){
  using namespace IMP::npctransport;
  using namespace IMP::atom;
  using namespace IMP::display;

  Hierarchy root( sd.get_root() );
  Hierarchies chains( get_fg_chains( root ) );
  unsigned int n_chains = chains.size();
  for(unsigned int i = 0 ; i < n_chains; i++) {
    Color color;
    // choose color
    if(n_chains <= 11) {
      color = get_display_color(i);
    } else {
      double f = i / (float)(n_chains - 1); // spread in [0..1]
      color = get_jet_color( f );
    }
    // apply color
    Hierarchies children = chains[i].get_children();
    for(unsigned int j = 0 ; j < children.size(); j++)
      {
        if( !Colored::particle_is_instance( children[i] ) ) {
          Colored::setup_particle( children[i], color );
        }
        else {
          Colored( children[i] ).set_color( color );
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
void set_fgs_in_cylinder( IMP::npctransport::SimulationData& sd, int n_layers )
{
  using namespace IMP;
  using atom::Hierarchy;
  using atom::Hierarchies;

  IMP::algebra::Cylinder3D cyl = sd.get_cylinder();
  Hierarchy root = sd.get_root() ;
  Hierarchies chains = IMP::npctransport::get_fg_chains(root);
  // compute the relative radius in which particles would be positioned
  // TODO: we assume here that particle radius is smaller
  //       than the cylinder radius - verify in runtime?
  double particle_radius =
    IMP::core::XYZR( chains[0].get_child(0) ).get_radius();
  // compute fraction of particle from full cylinder radius
  double relative_r =
    ( cyl.get_radius() - particle_radius ) / cyl.get_radius();
  // compute vertical poisition along central axis, and inter-layer distance
  double h_bottom_layer;
  if(n_layers == 1)
    h_bottom_layer = 0.5;
  else
    h_bottom_layer = 0.0;
  double dLayers = 0.0;
  if(n_layers > 1){
    dLayers = 1.0 / (n_layers - 1);
  }
  // calculate angle increments between adjacent fg nups in each layers
  int chains_per_layer =
    (int)( std::ceil( chains.size() / (n_layers + 0.0) ) );
  double angle_increments = 2 * IMP::PI / chains_per_layer;
  // pin chains to each layer
  for(int layer_num = 0 ; layer_num < n_layers ; layer_num++)
    {
      double relative_h =
        h_bottom_layer + layer_num * dLayers;
      for(int k = 0 ; k < chains_per_layer ; k++)
        {
          unsigned int chain_num = layer_num * chains_per_layer + k;
          if( chain_num >=  chains.size() )
            break; // may happen if len(chains) does not divide by n_layers
          double angle = k * angle_increments;
          algebra::Vector3D new_anchor =
            cyl.get_inner_point_at( relative_h, relative_r, angle);
          Hierarchy cur_chain( chains[chain_num] );
          core::XYZ d( cur_chain.get_child(0) );
          d.set_coordinates( new_anchor );
          d.set_coordinates_are_optimized( false);
          std::cout << "d = " << d << std::endl;
        }
    }
}


int main(int argc, char *argv[])
{
  // logging stuff:
  IMP::base::CreateLogContext clc("iteration");
  IMP::set_log_level(IMP::PROGRESS);
  RMF::set_show_hdf5_errors( true );
  // preparation:
  IMP_NPC_START_INT; // cmd-line options etc.
  IMP_NPC_PRINTHELP; // if help flag specified
  int num=IMP::npctransport::assign_ranges   // process config file
    (FLAGS_configuration, FLAGS_assignments,
     FLAGS_work_unit, FLAGS_show_steps);
  if(FLAGS_show_number_of_work_units)
    std::cout << "work units " << num;
  IMP::base::Pointer<IMP::npctransport::SimulationData> sd =
    new IMP::npctransport::SimulationData(FLAGS_assignments, FLAGS_statistics,
                                          FLAGS_quick, FLAGS_rmf_file);
  if(FLAGS_cylinder_anchoring)
    set_fgs_in_cylinder(*sd, 4);
  color_fgs( *sd );
  boost::timer timer= IMP::npctransport::create_boost_timer();
  int ntrials = sd->get_number_of_trials();
  // run:
  for(int i = 0 ; i < ntrials ; i++){
    std::cout <<  "Trial #" << i << std::endl;
    IMP::set_log_level(IMP::SILENT);
    if(! FLAGS_quick) // TODO: why is that?
      sd->reset_rmf();
    std::cout << "Initializing..." << std::endl;
    IMP::npctransport::initialize_positions(sd);
    std::cout << "Running..." << std::endl;
    sd->get_bd()->set_log_level(IMP::PROGRESS);
    IMP::benchmark::Profiler p;
    if(FLAGS_profile){
      p.set("profiling.pprof");
      std::cout << "Profiling begins..." << std::endl;
    }
    sd->get_bd()->optimize( sd->get_number_of_frames() );
    if(FLAGS_profile) {
      p.reset();
      std::cout << "Profiling end..." << std::endl;
    }
    sd->update_statistics( timer );
    std::cout << "Writing..." << std::endl;
    sd->write_geometry( FLAGS_final_configuration );
  }

  return 0;
 }
