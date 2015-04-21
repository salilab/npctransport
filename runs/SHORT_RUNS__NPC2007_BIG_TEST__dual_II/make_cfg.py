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
rest_length_factor = 1.25 #
obstacle_inflate_factor = 1.3
fg_coarse_factor=3.0 #
# fetch params from cmd-line
if(len(sys.argv)<=1):
    print " Bad Usage"
    exit(-1)
outfile = sys.argv[1]
if len(sys.argv) > 2:
    kaps_R = float(sys.argv[2])
print "kaps_R = %.2f" % (kaps_R)
if len(sys.argv) > 3:
    fg_coarse_factor=float(sys.argv[3])
print "coarse_factor = %.2f" % fg_coarse_factor
print "rest_length_factor = %.2f" % fg_coarse_factor

S=1

def get_basic_config():
    global S
    config = Configuration()
    IMP.npctransport.set_default_configuration(config)
    config.statistics_fraction.lower=1.0
    config.interaction_k.lower=10
    config.interaction_range.lower=1
    config.backbone_k.lower=0.25
    config.time_step_factor.lower=0.5 #### NOTE THIS ####
    ################## ACTIVE RANGE ####################
    create_range(config.time_step_factor, .5, 1.5, 3, 1)
    S=S*3
    config.time_step_wave_factor.lower=10.0 ### NOTE THIS ###
    ################## ACTIVE RANGE ####################
    create_range(config.time_step_wave_factor, 1, 10, 2, 1)
    S=S*2

    config.excluded_volume_k.lower=2
    config.nonspecific_range.lower=4
    config.nonspecific_k.lower=0.1
    ################## ACTIVE RANGE ####################
    create_range(config.nonspecific_k, .1, 1, 4, 1)
    S=S*4
    config.slack.lower = 7.5
    config.number_of_trials=1
    config.dump_interval_ns=2
    config.simulation_time_ns=5000
    config.angular_D_factor.lower=0.05 #lower to account for increased dynamic viscosity
                                      # in crowded environment and for coarse graining
    ################## ACTIVE RANGE ####################
    create_range(config.angular_D_factor, .05, 1, 4, 2)
    S=S*4
    config.statistics_interval_ns=0.1
    config.fg_anchor_inflate_factor=3.0/math.sqrt(fg_coarse_factor)
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



# ************** MAIN: *************
IMP.set_log_level(IMP.SILENT)
config= get_basic_config()

# Add floaters
n_kap_interactions=12/3
kaps= IMP.npctransport.add_float_type(config,
                                     number=3,
                                     radius=kaps_R,
                                      interactions= n_kap_interactions,
                                      type_name="kap")
############### ACTIVE RANGE #############
create_range(kaps.interaction_k_factor, lb=0.001, ub=10, steps = 10, base=1.5)
S=S*10
############### ACTIVE RANGE #############
create_range(kaps.interaction_range_factor, lb=0.5, ub=2, steps = 4, base=1)
S=S*4
############### ACTIVE RANGE #############
create_range(kaps.interactions, lb=4, ub=12, steps = 3, base=1.5)
S=S*3
#create_range(kaps.radius, lb = 10, ub = 30, steps = 5, base = 1)
nonspecifics= IMP.npctransport.add_float_type(config,
                                              number=3,
                                              radius=kaps_R, #-1,
                                              interactions=0,
                                              type_name="crap0")

fg_type="fg"
nbeads = 36
coarse_nbeads = 1 + int(math.ceil(nbeads / fg_coarse_factor)) # +1 for anchor
fg= IMP.npctransport.add_fg_type(config,
                                  type_name= fg_type,
                                  number_of_beads=coarse_nbeads,
                                  number=1,
                                  radius=7 * math.sqrt(fg_coarse_factor),
                                  interactions= int(math.ceil(1 * fg_coarse_factor)),
                                  rest_length_factor = rest_length_factor)
#    create_range(fg.number_of_beads,1,coarse_nbeads,steps=5,base=1)
##################  ACTIVE RANGE #############
create_range(fg.number,1,3,steps=3,base=1)
S=S*3
##################  ACTIVE RANGE #############
create_range(fg.rest_length_factor,1,1.5,steps=3,base=1)
S=S*3
add_interactions_for_fg(fg_type, k_fgkap)



# add bounding volumes
config.box_is_on.lower=1
config.box_side.lower=400
config.slab_is_on.lower=0
#config.tunnel_radius.lower=#
#config.slab_thickness.lower=250.0 # yeast nuclear envelope - see http://books.google.com/books?id=GvxdK1mdqQwC&pg=PA278&lpg=PA278&dq=yeast+nuclear+envelope+dimensions+nanometer&source=bl&ots=tHQoLfXHI1&sig=nRgZmLYnKuiRNP8n6vhm3bapjpI&hl=en&sa=X&ei=VtwKUtvAAsTAyAHOmIDYBg&ved=0CHsQ6AEwCA#v=onepage&q=yeast%20nuclear%20envelope%20dimensions%20nanometer&f=false
# config.slab_thickness.lower = max_z - config.fgs[0].radius.lower  # or also upper when there's steps?
#config.are_floaters_on_one_slab_side = 1 # all on top side


#create_range(interactionFG_FG.interaction_k, lb = 2.5, ub = 7.5, steps = 10, base = 1)
# internal FG-FG
#n_fg_types = len(config.fgs)
interactionFG_FG= IMP.npctransport.add_interaction(config,
                                                   name0= fg.type,
                                                   name1= fg.type,
                                                   interaction_k= k_fgfg,
                                                   interaction_range= range_fgfg)
############### ACTIVE RANGE #############
create_range(interactionFG_FG.interaction_k, lb=0.001, ub=10, steps = 10, base=1.33)
S=S*10
############### ACTIVE RANGE #############
create_range(interactionFG_FG.interaction_range, lb=0.5, ub=2, steps = 4, base=1)
S=S*4

print "S=",S
# dump to file
f=open(outfile, "wb")
f.write(config.SerializeToString())
# dump text
outfile_txt = os.path.splitext(outfile)[0] + ".txt"
f_txt = open(outfile_txt, "w")
print >>f_txt, config
f_txt.close()
