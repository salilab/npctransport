#!/usr/bin/python
from IMP.npctransport import *
import sys
import math
import os

# fetch params
# Usage: <cmd> <outfile> [kaps_R=25.0]
outfile = sys.argv[1]
kaps_R = 15.0
fg_coarse_factor=3 # 2
k_fgfg=1.0
k_fgkap=5.0
rest_length_factor = 1.2 # 1
obstacle_inflate_factor = 1.3

if(len(sys.argv) > 2):
    kaps_R = float(sys.argv[2])
print "kaps_R = %.2f" % (kaps_R)
if(len(sys.argv) > 3):
    fg_coarse_factor = float(sys.argv[3])
print "fg_coarse_factor = %.1f" % (fg_coarse_factor)

def get_basic_config():
    config = Configuration()
    IMP.npctransport.set_default_configuration(config)
    config.statistics_fraction.lower=1.0
    #config.dump_interval=1
    config.interaction_k.lower=1
    config.interaction_range.lower=1
    # create_range(config.backbone_k, .2, 1, 10
    config.backbone_k.lower=0.25
    #config.time_step_factor.lower=0.3
    config.time_step_factor.lower=1
    #create_range(config.rest_length_factor, .5, 1, 10)
    config.excluded_volume_k.lower=2
    config.nonspecific_range.lower=4
    config.nonspecific_k.lower=0.1
    config.slack.lower = 7.5
    config.number_of_trials=1
    config.dump_interval_ns=1
    config.simulation_time_ns=5000
    config.angular_D_factor.lower=0.1 #increased dynamic viscosity relative to
                                      # water?
    config.statistics_interval_ns=0.1
    config.fg_anchor_inflate_factor=3/fg_coarse_factor
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
                                     interaction_range=4)
    if(k_kap_steps > 1):
        create_range(interactionFG_KAP.interaction_k,
                     k_kap_lower, k_kap_upper,
                     steps = k_kap_steps,
                     base = k_kap_base)
    interactionFG_CRAP= IMP.npctransport.add_interaction(config,
                                       name0=fg_name,
                                       name1="inert",
                                       interaction_k=0,
                                       interaction_range=0)

def add_fg(config, k, nbeads, type_name,
                    coarse_factor=fg_coarse_factor):
    coarse_nbeads = 1 + int(math.ceil(nbeads / coarse_factor)) # +1 for anchor
    fgs= IMP.npctransport.add_fg_type(config,
                                      type_name=type_name,
                                      number_of_beads= coarse_nbeads,
                                      number=k,
                                      radius=6 * math.sqrt(coarse_factor),
                                      interactions= int(math.ceil(1 * coarse_factor)),
                                      rest_length_factor = rest_length_factor)
    add_interactions_for_fg(type_name, k_fgkap)




# ************** MAIN: *************
IMP.set_log_level(IMP.SILENT)
config= get_basic_config()

add_fg(config, k=2, nbeads = 30, type_name="fg")

# add bounding volumes
config.box_is_on.lower=1
config.box_side.lower=200

# Add floaters
n_kap_interactions=12/3
kaps= IMP.npctransport.add_float_type(config,
                                     number=10,
                                     radius=kaps_R,
                                      interactions= n_kap_interactions,
                                      type_name="kap")
#create_range(kaps.interaction_k_factor, lb=1, ub=5, steps = 10, base=1)
#create_range(kaps.radius, lb = 10, ub = 30, steps = 5, base = 1)
nonspecifics= IMP.npctransport.add_float_type(config,
                                              number=10,
                                              radius=kaps_R, #-1,
                                              interactions=0,
                                              type_name="inert")
#create_range(nonspecifics.radius, lb = 10, ub = 30, steps = 5, base = 1)
# fg with kaps / craps
#add_interactions_for_fg("fg0", 2.5, 7.5, k_kap_steps = 10, k_kap_base=1)
#add_interactions_for_fg("fg1", 2.5, 7.5, k_kap_steps = 10, k_kap_base=1)

# non-specific attraction
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
                                                           interaction_k= k_fgfg,
                                                           interaction_range= 4)

# dump to file
f=open(outfile, "wb")
f.write(config.SerializeToString())
f.close()
# dump text
outfile_txt = os.path.splitext(outfile)[0] + ".txt"
f_txt = open(outfile_txt, "w")
print >>f_txt, config
f_txt.close()
