#!/usr/bin/python
#/**
# * \file fg_array.cpp
# * \brief Simulate an fg and a kap interacting
# * Copyright 2007-2012 IMP Inventors. All rights reserved.
# */

#define IMP_NPC_MAIN
#include <IMP/npctransport/main.h>

import RMF
import IMP.algebra
import IMP.core
import IMP.atom
import IMP.display
import IMP
from IMP import container
import IMP.base
import IMP.npctransport
import IMP.benchmark
import IMP.example
from optparse import OptionParser
import math

def get_cmdline_options(args = None):
    usage = \
        "Usage: %prog [options]" + \
        "\n\nRuns FG Nups simulation with a cylinder"
    parser = OptionParser(usage)
    parser.add_option("-w", "--work_unit",
                      type="int", default=1,
                      help="The work unit")
    parser.add_option("-c", "--configuration", metavar="CONFIG_FILE",
                      type="string", default="configuration.pb",
                      help="input configuration file in protobuf format" +
                      " [default: %default]")
    parser.add_option("-a", "--assignments", metavar="ASSIGN_FILE",
                      type="string", default="assignments.pb",
                      help="output assignments file in protobuf format" +
                      ", recording the assignment being executed" +
                      " [default: %default]")
    parser.add_option("-s", "--statistics", metavar="STATISTICS_FILE",
                      type="string", default="statistics.pb",
                      help="output statistics file in protobuf format" +
                      " [default: %default]")
    parser.add_option("-f", "--final_configuration", metavar="FINAL_FILE",
                      type="string", default="final.pym",
                      help="output final configuration file" +
                      " [default: %default]")
    parser.add_option("-r", "--rmf_file", metavar="RMF_FILE",
                      type="string", default="output.rmf",
                      help="RMF file for recording the simulation progress" +
                      " [default: %default]")
    parser.add_option("-p", "--profile",
                      action="store_true", default=False,
                      help="Whether to turn on profiling for the first run")
    parser.add_option("-q", "--quick",
                      action="store_true", default=False,
                      help="Reduce all steps to the minimum");
    parser.add_option("--show_steps",
                      action="store_true", default=False,
                      help="Show the steps for each modified variable");
    parser.add_option("--show_number_of_work_units",
                      action="store_true", default=False,
                      help="Show the number of work units");
    parser.add_option("--cylinder_anchoring",
                      action="store_true", default=False,
                      help="anchor FG nups to a cylinder specified in" +
                      "the config file")
    (options, args) = parser.parse_args()
    return options


def set_fg_grid( sd ):
    """
    anchors all the fgs to a planar surface
    at the edge of the simulation bounding box

    sd - SimulationData object
    """
    # get bottom surface of the simulation data bounding box
    lower_corner_XY =  IMP.algebra.Vector2D(
        sd.get_box().get_corner(0)[0],
        sd.get_box().get_corner(0)[1])
    upper_corner_XY =  IMP.algebra.Vector2D(
        sd.get_box().get_corner(1)[0],
        sd.get_box().get_corner(1)[1])
    surface = IMP.algebra.BoundingBox2D(
        lower_corner_XY, upper_corner_XY )
    # get fg
    root = sd.get_root() # atom.Hierarchy
    chains = IMP.npctransport.get_fg_chains(root) # atom.Hierarchies
    # anchor fgs to surface,
    # for now using random non-overlapping points
    r= IMP.core.XYZR( chains[0].get_child(0)).get_radius()
    sites = IMP.algebra.Vector2Ds()
    while ( len(sites) < len(chains) ):
        cur = IMP.algebra.get_random_vector_in(surface) # Vector2D
        bad = False;
        for i in range( len(sites) ):
            # 2*r non-overlapping
            if (IMP.algebra.get_distance(sites[i], cur) < 2*r) :
                bad=True
                break
        if (not bad):
            sites.append(cur)
            print "Site # ",len(sites), " is ", cur
    for i in range( len(chains) ):
        r = IMP.atom.Hierarchy( chains[i] )
        d = IMP.core.XYZ( r.get_child(0) )
        d.set_coordinates(IMP.algebra.Vector3D(
                sites[i][0],
                sites[i][1],
                sd.get_box().get_corner(0)[2]))
        d.set_coordinates_are_optimized(False);
        print "d = ", d


def color_fgs( sd ):
    """
    color the different fgs in different colors

    chains - the SimulationData object
    """
    root = sd.get_root() # atom.Hierarchy
    chains = IMP.npctransport.get_fg_chains( root ) # atom.Hierarchies
    n_chains = len( chains )
    for i in range( n_chains ):
        color = None
        if(n_chains <= 11): # suitable for small number of colors
            color = IMP.display.get_display_color(i)
        else: # large number of colors
            f = i / float(n_chains - 1) # spread in [0..1]
            color = IMP.display.get_jet_color( f )
        for child in chains[i].get_children():
            if(not IMP.display.Colored.particle_is_instance( child ) ):
                IMP.display.Colored.setup_particle( child, color )
            else:
                IMP.display.Colored( child ).set_color( color )

def set_fgs_in_cylinder( sd, n_layers ):
    """
    anchors the FGs to the surface of the simulation bounding cylinder
    (= slab constraint)

    sd - the SimulationData object
    n_layers - number of fg nup layers
    """
    cyl = sd.get_cylinder()
    root = sd.get_root() # atom.Hierarchy
    chains = IMP.npctransport.get_fg_chains(root) # atom.Hierarchies
    # compute the relative radius in which particles would be positioned
    # TODO: we assume here that particle radius is smaller
    #       than the cylinder radius - verify in runtime?
    particle_radius = \
        IMP.core.XYZR( chains[0].get_child(0)).get_radius()
    # compute fraction of particle from full cylinder radius
    relative_r = \
        ( cyl.get_radius() - particle_radius ) / cyl.get_radius()
    # compute vertical poisition along central axis, and inter-layer distance
    bottom_layer_height = None
    if(n_layers == 1):
        bottom_layer_height = 0.5
    else:
        bottom_layer_height = 0.0
    delta_layers = 0.0
    if(n_layers > 1):
        delta_layers = 1.0 / (n_layers - 1)
    # calculate angle increments between adjacent fg nups in each layers
    chains_per_layer = int( math.ceil( len(chains) / float(n_layers) ) )
    angle_increments = 2 * math.pi / chains_per_layer
    # pin chains to each layer
    for layer in range( n_layers ):
        relative_h = bottom_layer_height + layer * delta_layers
        for k in range( chains_per_layer ):
            chain_num = layer * chains_per_layer + k #
            if( chain_num >=  len(chains) ):
                break; # may happen if len(chains) does not divide by n_layers
            angle = k * angle_increments
            new_anchor = cyl.get_inner_point_at( \
                relative_h, relative_r, angle)
            cur_chain = IMP.atom.Hierarchy( chains[chain_num] )
            d = IMP.core.XYZ( cur_chain.get_child(0) )
            d.set_coordinates( new_anchor )
            d.set_coordinates_are_optimized(False);
            print "d = ", d


################## MAIN ####################
flags = get_cmdline_options()
print flags
# process info from protobuf
IMP.npctransport.assign_ranges(flags.configuration, flags.assignments,
                               flags.work_unit, flags.show_steps)
print "here"
if(flags.show_number_of_work_units):
    print "work units ", flags.work_units
RMF.set_show_hdf5_errors(True)
print "here2"
# #ifdef IMP_BENCHMARK_USE_GOOGLE_PERFTOOLS_PROFILE
# #define IMP_NPC_SET_PROF(tf) if (FLAGS_profile && i==0) {               \
#     IMP::benchmark::set_is_profiling(tf);                               \
#   }
# #else
# #define IMP_NPC_SET_PROF(tf)
# #endif

timer = IMP.npctransport.create_boost_timer()
print "here3"
IMP.set_log_level(IMP.WARNING)
sd = IMP.npctransport.SimulationData(
    flags.assignments,
    flags.statistics,
    flags.quick,
    flags.rmf_file)
print "here4"
print "RMF file: ", sd.get_rmf_file_name()
if(flags.cylinder_anchoring):
    set_fgs_in_cylinder(sd, 4)
color_fgs( sd )
ntrials = sd.get_number_of_trials()
print "Number of trials: ", ntrials
for i in range(ntrials):
    clc = IMP.base.CreateLogContext("iteration")
    if(not flags.quick): # TODO: why is that?
        sd.reset_rmf()
    print "Initializing..."
    IMP.npctransport.initialize_positions(sd)
    print "Running..."
    sd.get_bd().set_log_level(IMP.WARNING)
    p = IMP.benchmark.Profiler()
    if(flags.profile):
        p.start("profile.pprof")
        print "Profiling begins..."
    sd.get_bd().optimize( sd.get_number_of_frames() )
    if(flags.profile):
        p.stop()
        print "Profiling ends"
    sd.update_statistics( timer ) # TODO: timer?
    print "Writing..."
    sd.write_geometry( flags.final_configuration )
    # Profiling?
