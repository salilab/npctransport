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
#import IMP.benchmark
#import IMP.example
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
    parser.add_option("-i", "--initialization_rmf_file",
                      help="initialize the simulation from the last frame of " +
                      "the specified RMF file, assuming it was created by " +
                      "this simualtion process")
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


def get_fgs_of_type( fg_type, root ):
    """
    returns all fg chain hierarchies of type fg_type_name that descend from root

    fg_type - ParticleType object with the required type of fgs
    """
    ret = IMP.atom.Hierarchies()
    # I. return root itself if the type of its first direct child
    # is fg_type_name
    if (root.get_number_of_children() >0):
        c= root.get_child(0) # atom.Hierarchy
        if (IMP.core.Typed.particle_is_instance(c)):
            t=IMP.core.Typed(c).get_type()
            if (t == fg_type):
                ret.append( root )
    # II. Otherwise, recurse on all of root's children
    for i in range( root.get_number_of_children()) :
        ichild_fgs = get_fgs_of_type( fg_type, root.get_child(i) )
        ret.extend( ichild_fgs )
    return ret


# TODO: implement
def set_specific_fgs_in_cylinder( sd, fgs_list, n_layers,
                                  relative_bottom, relative_top):
    """
    anchors the FG chains in the list fgs_list to the surface of
    the simulation bounding cylinder (= tunnel inside slab).
    The chains are distributed as evenly as possible in n_layers,
    between relative_bottom and relative_top height (values from 0 to 1,
    which correspond to slab bottom and slab top respectively)

    sd - the SimulationData object
    fgs_list - list of objects of type Hierarchy, each supposed to be an fg
    n_layers - number of fg nup layers
    relative bottom - bottom layer position relative to the cylinder z-axis
                      (0 = cylinder bottom ; 1 = cylinder top)
    relative top - top layer position relative to the cylinder z-axis
                      (0 = cylinder bottom ; 1 = cylinder top)
    """
    cyl = sd.get_cylinder()
    # compute the relative radius in which particles would be positioned
    # TODO: we assume here that particle radius is smaller
    #       than the cylinder radius - verify in runtime?
    particle_radius = \
        IMP.core.XYZR( fgs_list[0].get_child(0)).get_radius()
    # compute fraction of particle from full cylinder radius
    relative_r = \
        ( cyl.get_radius() - particle_radius ) / cyl.get_radius()
    # compute vertical poisition along central axis, and inter-layer distance
    bottom_layer_height = None
    relative_mid = (relative_bottom + relative_top) / 2
    if(n_layers == 1):
        bottom_layer_height = relative_mid
    else:
        bottom_layer_height = relative_bottom
    delta_layers = 0.0
    if(n_layers > 1):
        delta_layers = (relative_top - relative_bottom ) / (n_layers - 1.0)
        print "Delta layers: %f (relative value)" % delta_layers
    # calculate angle increments between adjacent fg nups in each layers
    chains_per_layer = int( math.ceil( len(fgs_list) / float(n_layers) ) )
    angle_increments = 2.0 * math.pi / chains_per_layer
    # pin chains to each layer
    for layer in range( n_layers ):
        relative_h = bottom_layer_height + layer * delta_layers
        angle_phase = layer * angle_increments / n_layers
        for k in range( chains_per_layer ):
            chain_num = layer * chains_per_layer + k #
            if( chain_num >=  len(fgs_list) ):
                break; # may happen if len(chains) does not divide by n_layers
            angle = k * angle_increments + angle_phase
            new_anchor = cyl.get_inner_point_at( \
                relative_h, relative_r, angle)
            cur_chain = IMP.atom.Hierarchy( fgs_list[chain_num] )
            d = IMP.core.XYZ( cur_chain.get_child(0) )
            d.set_coordinates( new_anchor )
            d.set_coordinates_are_optimized(False);
            print "d = ", d

def set_fgs_in_cylinder( sd, n_layers ):
    """
    anchors the FGs to the surface of the simulation bounding cylinder
    (= slab constraint)

    sd - the SimulationData object
    n_layers - number of fg nup layers
    """
    cyl = sd.get_cylinder()
    root = sd.get_root() # atom.Hierarchy
    fg_chains = IMP.npctransport.get_fg_chains(root) # atom.Hierarchies
    set_specific_fgs_in_cylinder(sd = sd,
                        fgs_list = fg_chains,
                        n_layers = n_layers,
                        relative_bottom = 0.0, relative_top = 1.0)

def set_fgs_three_types( sd ):
    """
    anchors the FGs to the surface of the simulation bounding cylinder
    (= slab constraint), using fg0 for top filaments, fg1 for middle,
    and fg2 for bottom

    sd - the SimulationData object
    """
    cyl = sd.get_cylinder()
    root = sd.get_root() # atom.Hierarchy
    fgs_cyt = get_fgs_of_type(IMP.npctransport.get_type_of_fg(0), root)
    fgs_middle = get_fgs_of_type(IMP.npctransport.get_type_of_fg(1), root)
    fgs_nuclear = get_fgs_of_type(IMP.npctransport.get_type_of_fg(2), root)
    print "DEBUG stats:"
    print fgs_cyt
    print fgs_cyt[0]
    print type(fgs_cyt[0])

    set_specific_fgs_in_cylinder(sd = sd,
                        fgs_list = fgs_cyt,
                        n_layers = 1,
                        relative_bottom = 1.0, relative_top = 1.0)
    set_specific_fgs_in_cylinder(sd = sd,
                        fgs_list = fgs_middle,
                        n_layers = 3,
                        relative_bottom = 0.2, relative_top = 0.8)
    set_specific_fgs_in_cylinder(sd = sd,
                        fgs_list = fgs_nuclear,
                        n_layers = 1,
                        relative_bottom = 0.0, relative_top = 0.0)

def optimize_in_chunks( sd, nframes, nchunks ):
    """
    Optimizes sd->bd() in nchunks iterations, writing statistics at
    the end of each iteration
    """
    timer = IMP.npctransport.create_boost_timer()
    nframes_left = nframes
    nframes_chunk= math.ceil(nframes_left / nchunks)
    while(nframes_left > 0):
        nframes_chunk = min(nframes_chunk, nframes_left)
        sd.get_bd().optimize( nframes_chunk )
        sd.update_statistics( timer ) # TODO: timer?
        nframes_left = nframes_left - nframes_chunk


################## MAIN ####################
flags = get_cmdline_options()
print flags
# process info from protobuf, using [work_unit]'th combination of values
n = IMP.npctransport.assign_ranges(flags.configuration, flags.assignments,
                               flags.work_unit, flags.show_steps)
if(flags.show_number_of_work_units):
    print "total number of work units ", n
RMF.set_show_hdf5_errors(True)
# #ifdef IMP_BENCHMARK_USE_GOOGLE_PERFTOOLS_PROFILE
# #define IMP_NPC_SET_PROF(tf) if (FLAGS_profile && i==0) {               \
#     IMP::benchmark::set_is_profiling(tf);                               \
#   }
# #else
# #define IMP_NPC_SET_PROF(tf)
# #endif

IMP.set_log_level(IMP.PROGRESS)
sd = IMP.npctransport.SimulationData(
    flags.assignments,
    flags.statistics,
    flags.quick,
    flags.rmf_file)
print "RMF file: ", sd.get_rmf_file_name()
print get_fgs_of_type(IMP.npctransport.get_type_of_fg(0), sd.get_root())
if(flags.cylinder_anchoring):
#    set_fgs_in_cylinder(sd, 4)
    set_fgs_three_types(sd)
color_fgs( sd )
ntrials = sd.get_number_of_trials()
print "Number of trials: ", ntrials
for i in range(ntrials):
    clc = IMP.base.CreateLogContext("iteration")
    if(not flags.quick): # TODO: why is that?
        sd.reset_rmf()
    print "Initializing..."
    if(flags.initialization_rmf_file):
        sd.initialize_positions_from_rmf(
            flags.initialization_rmf_file )
    else:
        IMP.npctransport.initialize_positions( sd )
    nframes_running = math.ceil(sd.get_statistics_fraction()
                                * sd.get_number_of_frames())
    nframes_equilib = sd.get_number_of_frames() - nframes_running
    print "Equilibrating for frames...", nframes_equilib
    sd.reset_rmf() # TODO: should I keep it? it means no RMF for initialization
    sd.get_bd().optimize(nframes_equilib)
    print "Running for frames...", nframes_running
    sd.reset_statistics_optimizer_states(); # only take stats after equib
    sd.get_bd().set_log_level(IMP.PROGRESS)
    p = IMP.benchmark.Profiler()
    if(flags.profile):
        p.start("profile.pprof")
        print "Profiling begins..."
    sd.get_bd().set_current_time(0)
    nchunks= 2500 # parametrize externally?
    optimize_in_chunks(sd, nframes_running, nchunks)
    if(flags.profile):
        p.stop()
        print "Profiling ends"
    print "Writing..."
    sd.write_geometry( flags.final_configuration )
    # Profiling?
