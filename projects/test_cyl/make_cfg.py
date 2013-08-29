#!/usr/bin/python
from IMP.npctransport import *
import sys
import math
import re
import os

# defaults
kaps_R = 30.0
k_fgfg=0.25
range_fgfg=4.0
k_fgkap=3.0
range_fgkap=4.0
rest_length_factor = 1.25 # 1
obstacle_inflate_factor = 1.3
fg_coarse_factor=3.0 # 3
# fetch params from cmd-line
if(len(sys.argv)<=1):
    print " Usage: <cmd> <outfile> [kaps_R=%.1f] [k_fgfg=%.1f]" % (kaps_R, k_fgfg)
    exit(-1)
outfile = sys.argv[1]
if len(sys.argv) > 2:
    kaps_R = float(sys.argv[2])
print "kaps_R = %.2f" % (kaps_R)
if len(sys.argv) > 3:
    k_fgfg = float(sys.argv[3])
print "k_fgfg = %.2f" %k_fgfg

def get_basic_config():
    config = Configuration()
    IMP.npctransport.set_default_configuration(config)
    config.statistics_fraction.lower=1.0
    config.interaction_k.lower=10
    config.interaction_range.lower=1
    config.backbone_k.lower=0.25
    config.time_step_factor.lower=1 #### NOTE THIS ####
    #create_range(config.rest_length_factor, .5, 1, 10)
    config.time_step_wave_factor.lower=5 #### NOTE THIS ####
    config.excluded_volume_k.lower=2
    config.nonspecific_range.lower=4
    config.nonspecific_k.lower=0.1
    config.slack.lower = 7.5
    config.number_of_trials=1
    config.dump_interval_ns=5
    config.simulation_time_ns=5000
    config.angular_D_factor.lower=0.05 #lower to account for increased dynamic viscosity
                                      # in crowded environment and for coarse graining
    config.statistics_interval_ns=0.1
    config.fg_anchor_inflate_factor=3.0/math.sqrt(fg_coarse_factor)
    config.is_exclude_floaters_from_slab_initially=0
    return config



def add_interactions_for_fg(fg_name,
                            k_kap_lower,
                            k_kap_upper = 0, # relevant only if k_kap_steps > 1
                            k_kap_steps = 1,
                            k_kap_base = 1):
    interactionFG_KAP= IMP.npctransport.add_interaction(config,
                                                        name0=fg_name,
                                                        name1="kap",
                                                        interaction_k=k_kap_lower,
                                                        interaction_range=range_fgkap)
    if(k_kap_steps > 1):
        create_range(interactionFG_KAP.interaction_k,
                     k_kap_lower, k_kap_upper,
                     steps = k_kap_steps,
                     base = k_kap_base)
    interactionFG_CRAP= IMP.npctransport.add_interaction(config,
                                       name0=fg_name,
                                       name1="crap0",
                                       interaction_k=0,
                                       interaction_range=0)

def add_fg_based_on(config, type_name, k, nbeads, origin=None,
                    coarse_factor=fg_coarse_factor):
    ''' Read mrc_filename, cluster to k clusters, and create k
        fgs with nbeads/coarse_factor beads, anchored at the clusters,
        normalized by mean_loc. An additional bead is used for the
        anchor residue

        @param coarse_factor - a factor by which to coarse grain beads
                               (e.g. 3 beads : 1)

        @note if origin==None, initiate it as a 3D coordinate that
        is the mean coordinate of the MRC file in mrc_filename.
        @note also updates global variables max_xy, max_x, max_y, max_z

        @return the mean location of the MRC file
        '''
    global max_r, max_x, max_y, max_z
    coarse_nbeads = 1 + int(math.ceil(nbeads / coarse_factor)) # +1 for anchor
    fgs= IMP.npctransport.add_fg_type(config,
                                      type_name= type_name,
                                      number_of_beads= coarse_nbeads,
                                      number=k,
                                      radius=7 * math.sqrt(coarse_factor),
                                      interactions= int(math.ceil(1 * coarse_factor)),
                                      rest_length_factor = rest_length_factor)
    add_interactions_for_fg(type_name, k_fgkap)

def add_obstacle(config, mrc_filename, k, R, origin=None):
    ''' Read mrc_filename, cluster to k clusters, and create k
        obstacles anchored at the clusters, normalized about origin.

        @param k the number of obstacles
        @param R the obstacle radius, to be inflated by obstacle_inflate_factor
        @param origin the new origin

        @note if origin==None, initiate it as a 3D coordinate that
        is the mean coordinate of the MRC file in mrc_filename.

        @return the mean location of the MRC file clusters
        '''
    # get type name as filename without folder and extension parts
    print mrc_filename
    type_search= re.search("([^/]*?)(?:[.].*)*$", mrc_filename)
    type_name = type_search.groups(0)[0]
    # cluster anchors from MRC file
    centers, mean_loc = read_nups.cluster_MRC_file_with_cache(mrc_filename, k)
    if(origin is None):
        origin = mean_loc
    obstacle = IMP.npctransport.add_obstacle_type \
        (config, type_name=type_name, R = R * obstacle_inflate_factor)
    for center in centers:
        pos=obstacle.xyzs.add()
        pos.x=center[0] - origin[0]
        pos.y=center[1] - origin[1]
        pos.z=center[2] - origin[2]
        r = math.sqrt(pos.x**2 + pos.y**2)
        print "OBSTACLE", mrc_filename, "z=", pos.z, "r=", r, "R=", R
    return mean_loc



# ************** MAIN: *************
IMP.set_log_level(IMP.base.SILENT)
config= get_basic_config()

# Add floaters
n_kap_interactions=12/3
kaps= IMP.npctransport.add_float_type(config,
                                     number=11,
                                     radius=kaps_R,
                                      interactions= n_kap_interactions,
                                      type_name="kap")
kaps.k_z_bias.lower=0.005
kaps.k_z_bias_fraction.lower=0.3
############### ACTIVE RANGE #############
create_range(kaps.interaction_k_factor, lb=1, ub=5, steps = 10, base=1)
##########################################
#create_range(kaps.radius, lb = 10, ub = 30, steps = 5, base = 1)
nonspecifics= IMP.npctransport.add_float_type(config,
                                              number=2,
                                              radius=kaps_R, #-1,
                                              interactions=0,
                                              type_name="crap0")

#create_range(nonspecifics.radius, lb = 10, ub = 30, steps = 5, base = 1)
# fg with kaps / craps
#add_interactions_for_fg("fg0", 2.5, 7.5, k_kap_steps = 10, k_kap_base=1)
#add_interactions_for_fg("fg1", 2.5, 7.5, k_kap_steps = 10, k_kap_base=1)

# non-specific attraction
#create_range(config.nonspecific_k,lb=1.0,ub=2.0,steps = 3,base=1)


# Add FGs with anchors
# (Stoicheometries from Alber et al. 2007b, Determining..., Fig. 3)
add_fg_based_on(config, "fg", k=1, nbeads=3)
# add bounding volumes
config.box_is_on.lower=1
config.box_side.lower=300
config.slab_is_on.lower=1
config.tunnel_radius.lower=150
config.slab_thickness.lower=100.0 # yeast nuclear envelope - see http://books.google.com/books?id=GvxdK1mdqQwC&pg=PA278&lpg=PA278&dq=yeast+nuclear+envelope+dimensions+nanometer&source=bl&ots=tHQoLfXHI1&sig=nRgZmLYnKuiRNP8n6vhm3bapjpI&hl=en&sa=X&ei=VtwKUtvAAsTAyAHOmIDYBg&ved=0CHsQ6AEwCA#v=onepage&q=yeast%20nuclear%20envelope%20dimensions%20nanometer&f=false
# config.slab_thickness.lower = max_z - config.fgs[0].radius.lower  # or also upper when there's steps?
config.are_floaters_on_one_slab_side = 1 # all on top side


#create_range(interactionFG_FG.interaction_k, lb = 2.5, ub = 7.5, steps = 10, base = 1)
# internal FG-FG
#n_fg_types = len(config.fgs)
for i,fg0 in enumerate(config.fgs):
    for j,fg1 in enumerate(config.fgs):
        if(i>j): continue # avoid duplicates
        interactionFG_FG= IMP.npctransport.add_interaction(config,
                                                           name0= fg0.type,
                                                           name1= fg1.type,
                                                           interaction_k= k_fgfg,
                                                           interaction_range= range_fgfg)

# dump to file
f=open(outfile, "wb")
f.write(config.SerializeToString())
# dump text
outfile_txt = os.path.splitext(outfile)[0] + ".txt"
f_txt = open(outfile_txt, "w")
print >>f_txt, config
f_txt.close()
