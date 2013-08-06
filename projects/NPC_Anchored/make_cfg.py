#!/usr/bin/python
from IMP.npctransport import *
import sys
import math
import read_nups
import re

# fetch params
# Usage: <cmd> <outfile> [kaps_R=25.0]
outfile = sys.argv[1]
kaps_R = 15.0
if(len(sys.argv) > 2):
    kaps_R = float(sys.argv[2])
print "kaps_R = %.2f" % (kaps_R)

def get_basic_config():
    config = Configuration()
    IMP.npctransport.set_default_configuration(config)
    config.statistics_fraction.lower=1.0
    #config.dump_interval=1
    config.interaction_k.lower=10
    config.interaction_range.lower=1
    # create_range(config.backbone_k, .2, 1, 10
    config.backbone_k.lower=1
    #config.time_step_factor.lower=0.3
    config.time_step_factor.lower=3
    #create_range(config.rest_length_factor, .5, 1, 10)
    config.excluded_volume_k.lower=20
    config.nonspecific_range.lower=2
    config.nonspecific_k.lower=0.01
    config.slack.lower = 15
    config.number_of_trials=1
    config.dump_interval_ns=0.1
    config.simulation_time_ns=1000
    config.angular_D_factor.lower=0.3 #increased dynamic viscosity relative to water?
    config.statistics_interval_ns=0.1
    ###
    #simulation bounding volumes:
    config.box_is_on.lower=1
    config.box_side.lower=400
    config.slab_is_on.lower=0
    config.slab_thickness.lower=150
    config.tunnel_radius.lower=75
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
                                     interaction_range=2)
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

def add_fg_based_on(config, mrc_filename, k, nbeads, mean_loc=None):
    ''' Read mrc_filename, cluster to k clusters, and create k
        fgs with nbeads beads, anchored at the clusters, normalized by
        mean_loc.

        @note if mean_loc==None, initiate it as a 3D coordinate that
        is the mean coordinate of the MRC file in mrc_filename.
        @note also updates global variables max_xy, max_x, max_y, max_z

        @return the mean location of the MRC file
        '''
    global max_r, max_x, max_y, max_z
    # get type name as filename without folder and extension parts
    type_search= re.search("([^/]*?)(?:[.].*)*$", mrc_filename)
    type_name = type_search.groups(0)[0]
    # cluster anchors from MRC file
    kmeans, centers, new_mean, pos_voxels = read_nups.cluster_MRC_file(mrc_filename, k)
    if(mean_loc is None):
        mean_loc = new_mean
    fgs= IMP.npctransport.add_fg_type(config,
                                      type_name= type_name,
                                      number_of_beads= nbeads,
                                      number=len(centers),
                                      radius=6,
                                      interactions=1,
                                      rest_length_factor = 1.5)
    add_interactions_for_fg(type_name, 3.0)
    for center in centers:
        pos=fgs.anchor_coordinates.add()
        pos.x=center[0] - mean_loc[0]
        pos.y=center[1] - mean_loc[1]
        pos.z=center[2] - mean_loc[2]
        print pos
        max_r = max(max_r, math.sqrt(pos.x**2 + pos.y**2) )
        max_x = max(max_x, abs(pos.x))
        max_y = max(max_y, abs(pos.y))
        max_z = max(max_z, abs(pos.z))
    return new_mean

def add_obstacle(config, mrc_filename, k, R, mean_loc=None):
    ''' Read mrc_filename, cluster to k clusters, and create k
        obstacles of radius R, anchored at the clusters, normalized by
        mean_loc.

        @note if mean_loc==None, initiate it as a 3D coordinate that
        is the mean coordinate of the MRC file in mrc_filename.

        @return the mean location of the MRC file
        '''
    # get type name as filename without folder and extension parts
    print mrc_filename
    type_search= re.search("([^/]*?)(?:[.].*)*$", mrc_filename)
    type_name = type_search.groups(0)[0]
    # cluster anchors from MRC file
    kmeans, centers, new_mean, pos_voxels = read_nups.cluster_MRC_file(mrc_filename, k)
    if(mean_loc is None):
        mean_loc = new_mean
    obstacle = IMP.npctransport.add_obstacle_type \
        (config, type_name=type_name, R=R)
    for center in centers:
        pos=obstacle.xyzs.add()
        pos.x=center[0] - mean_loc[0]
        pos.y=center[1] - mean_loc[1]
        pos.z=center[2] - mean_loc[2]
        print pos
    return new_mean



# ************** MAIN: *************
config= get_basic_config()
config.dump_interval_ns=0.1
config.simulation_time_ns=10

# Add FGs with anchors
# (Stoicheometries from Alber et al. 2007b, Determining..., Fig. 3)
max_r=0
max_x=0
max_y=0
max_z=0
mean_loc=(add_fg_based_on(config, "MRCs/Nup57_16copies_chimera.mrc", k=16, nbeads=16))
add_fg_based_on(config, "MRCs/Nup49_16copies.mrc", k=16, nbeads = 17, mean_loc=mean_loc)
add_fg_based_on(config, "MRCs/Nsp1_16copies_1.mrc", k=16, nbeads = 33, mean_loc=mean_loc)
add_fg_based_on(config, "MRCs/Nsp1_16copies_2.mrc", k=16, nbeads = 33, mean_loc=mean_loc)
add_fg_based_on(config, "MRCs/Nup159_8copies.mrc", k=8, nbeads=20, mean_loc=mean_loc) # nbeads 20-24 = real number for Nup159, depending how you count double motifs
add_fg_based_on(config, "MRCs/Nup116_8copies_chimera.mrc", k=8, nbeads=46, mean_loc=mean_loc)
add_fg_based_on(config, "MRCs/Nup42_8copies_chimera.mrc", k=8, nbeads=21, mean_loc=mean_loc) # nbeads 21-27, depending on treratment of double motifs
add_fg_based_on(config, "MRCs/Nup100_8copies_chimera.mrc", k=8, nbeads=44, mean_loc=mean_loc)
add_fg_based_on(config, "MRCs/Nup145N_8copies_1_chimera.mrc", k=8, nbeads=44, mean_loc=mean_loc)
add_fg_based_on(config, "MRCs/Nup145N_8copies_2_chimera.mrc", k=8, nbeads=44, mean_loc=mean_loc)
# Nuclear:
add_fg_based_on(config, "MRCs/Nup1_8copies.mrc", k=8, nbeads=27, mean_loc=mean_loc)
add_fg_based_on(config, "MRCs/Nup59_16copies.mrc", k=16, nbeads=5, mean_loc=mean_loc) # contains also RRM Nup35-type domain for RNA binding (res 265-394 ; Uniprot), which supposedly overlaps some of the FGs
add_fg_based_on(config, "MRCs/Nup60_8copies.mrc", k=8, nbeads=11, mean_loc=mean_loc) # nup60 is supposed to tether Nup2 (depending on Gsp1p-GTP (Ran) switch, for which Nup2 has a binding site 583-720 ; Denning, Rexach et al. JCB 2001) ; Nup2 also interacts with cargo 35-50 and RNA (RRM Nup35-type domain) - from Uniprot)
add_fg_based_on(config, "MRCs/Nup53_16copies_chimera.mrc", k=16, nbeads=4, mean_loc=mean_loc) # contains also RRM Nup35-type domain for RNA binding (res 247-352 ; Uniprot), which supposedly overlaps some of the FGs ; and also a PSE1/Kap121 binding domain in a non-FG fashion, for which there might be a crystal structure (405-438 ; Uniprot)

# Add Structural nups as obstacles
# (Alber et al. 2007b, Deteriming..., Figure 3)
add_obstacle(config, "MRCs/Nup192_16copies.mrc", k=16, R=40, mean_loc=mean_loc)
add_obstacle(config, "MRCs/Nup188_16copies.mrc", k=16, R=40, mean_loc=mean_loc)
add_obstacle(config, "MRCs/Nup170_16copies.mrc", k=16, R=40, mean_loc=mean_loc)
add_obstacle(config, "MRCs/Nup170_16copies.mrc", k=16, R=40, mean_loc=mean_loc)
add_obstacle(config, "MRCs/Nup157_16copies.mrc", k=16, R=40, mean_loc=mean_loc)
add_obstacle(config, "MRCs/Nup133_16copies.mrc", k=16, R=30, mean_loc=mean_loc)
add_obstacle(config, "MRCs/Nup120_16copies.mrc", k=16, R=30, mean_loc=mean_loc)
add_obstacle(config, "MRCs/Nic96_16copies_1.mrc", k=16, R=30, mean_loc=mean_loc)
add_obstacle(config, "MRCs/Nic96_16copies_2.mrc", k=16, R=30, mean_loc=mean_loc)
add_obstacle(config, "MRCs/Nup85_16copies.mrc", k=16, R=30, mean_loc=mean_loc)
add_obstacle(config, "MRCs/Nup84_16copies.mrc", k=16, R=30, mean_loc=mean_loc)
add_obstacle(config, "MRCs/Nup82_8copies_1.mrc", k=8, R=35, mean_loc=mean_loc)
add_obstacle(config, "MRCs/Nup82_8copies_2.mrc", k=8, R=35, mean_loc=mean_loc)
add_obstacle(config, "MRCs/Nup145C_16copies.mrc", k=16, R=30, mean_loc=mean_loc)
add_obstacle(config, "MRCs/Ndc1_16copies.mrc", k=16, R=30, mean_loc=mean_loc)
add_obstacle(config, "MRCs/Gle1_8copies.mrc", k=8, R=30, mean_loc=mean_loc)
add_obstacle(config, "MRCs/Gle2_16copies.mrc", k=16, R=23, mean_loc=mean_loc)
add_obstacle(config, "MRCs/Seh1_16copies.mrc", k=16, R=22, mean_loc=mean_loc)
add_obstacle(config, "MRCs/Pom34_16copies.mrc", k=16, R=20, mean_loc=mean_loc)
add_obstacle(config, "MRCs/Sec13_16copies.mrc", k=16, R=21, mean_loc=mean_loc)

# add bounding volumes
config.box_is_on.lower=1
config.box_side.lower=max(max_z,max_x,max_y)*6 # 2000
config.slab_is_on.lower=1
config.tunnel_radius.lower=max_r - config.fgs[0].radius.lower # or also upper when there's steps?
config.slab_thickness.lower=max_z - config.fgs[0].radius.lower  # or also upper when there's steps?

# Add floaters
kaps= IMP.npctransport.add_float_type(config,
                                     number=30,
                                     radius=kaps_R,
                                     interactions=12)
#create_range(kaps.radius, lb = 10, ub = 30, steps = 5, base = 1)
nonspecifics= IMP.npctransport.add_float_type(config,
                                              number=30,
                                              radius=kaps_R, #-1,
                                              interactions=0)
#create_range(nonspecifics.radius, lb = 10, ub = 30, steps = 5, base = 1)
# fg with kaps / craps
#add_interactions_for_fg("fg0", 2.5, 7.5, k_kap_steps = 10, k_kap_base=1)
#add_interactions_for_fg("fg1", 2.5, 7.5, k_kap_steps = 10, k_kap_base=1)

# non-specific attraction
config.nonspecific_range.lower= 1.5
#create_range(config.nonspecific_k,lb=1.0,ub=2.0,steps = 3,base=1)

#create_range(interactionFG_FG.interaction_k, lb = 2.5, ub = 7.5, steps = 10, base = 1)
# internal FG-FG
#n_fg_types = len(config.fgs)
for i,fg0 in enumerate(config.fgs):
    for j,fg1 in enumerate(config.fgs):
        if(i>j): continue # avoid duplicates
        interactionFG_FG= IMP.npctransport.add_interaction(config,
                                                           name0= fg0.type,
                                                           name1= fg1.type,
                                                           interaction_k= 3.5,
                                                           interaction_range= 2)

# dump to file
f=open(outfile, "wb")
f.write(config.SerializeToString())
print config
