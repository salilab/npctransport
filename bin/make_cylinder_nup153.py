#!/usr/bin/env python

from IMP.npctransport import *
import sys
import barak_basic_configuration

def add_interactions_for_fg(fg_name,
                            k_kap_lower,
                            k_kap_upper = 0, # relevant only if k_kap_steps > 1
                            k_kap_steps = 1):
    interactionFG_KAP= IMP.npctransport.add_interaction(config,
                                     name0=fg_name,
                                     name1="kap",
                                     interaction_k=k_kap_lower,
                                     interaction_range=2)
    if(k_kap_steps > 1):
        create_range(interactionFG_KAP.interaction_k,
                     k_kap_lower, k_kap_upper,
                     steps = k_kap_steps)
    interactionFG_CRAP= IMP.npctransport.add_interaction(config,
                                       name0=fg_name,
                                       name1="crap0",
                                       interaction_k=0,
                                       interaction_range=0)

# ********* MAIN: *********
config= barak_basic_configuration.get_basic_config()
config.dump_interval_ns=25
config.simulation_time_ns=10000
config.box_is_on.lower=1
config.box_side.lower=300
config.slab_is_on.lower=1
config.slab_thickness.lower=150
config.tunnel_radius.lower=90

fg_cyto= IMP.npctransport.add_fg_type(config,
                                 number_of_beads=8,
                                 number=5,
                                 radius=8,
                                 interactions=1)
fg_middle= IMP.npctransport.add_fg_type(config,
                                 number_of_beads=8,
                                 number=15,
                                 radius=8,
                                 interactions=1)
fg_nuclear= IMP.npctransport.add_fg_type(config,
                                 number_of_beads=8,
                                 number=5,
                                 radius=8,
                                 interactions=1)
kap= IMP.npctransport.add_float_type(config,
                                     number=10,
                                     radius=15,
                                     interactions=12)
nonspecifics= IMP.npctransport.add_float_type(config,
                                              number=20,
                                              radius=15,
                                              interactions=0)

# fg with kaps / craps
add_interactions_for_fg("fg0",5)
add_interactions_for_fg("fg1",
                        k_kap_lower = 2.5,
                        k_kap_upper = 20,
                        k_kap_steps = 10)
add_interactions_for_fg("fg2",5)
# non-specific attraction
create_range(config.nonspecific_k, 0.1, 10, steps=5)
create_range(config.nonspecific_range, 0.1, 2, steps=3)
# internal FG-FG
for i in range(3):
    for j in range(i,3):
        interactionFG_FG= IMP.npctransport.add_interaction(config,
                                                           name0="fg%d" % i,
                                                           name1="fg%d" % j,
                                                           interaction_k=1,
                                                           interaction_range=2)

# dump to file
f=open(sys.argv[1], "wb")
f.write(config.SerializeToString())
print config

# dump quick file version
if len(sys.argv)>2:
    config.number_of_trials.lower=2
    config.number_of_frames.lower=2
    config.dump_interval=1
    f=open(sys.argv[2], "wb")
    f.write(config.SerializeToString())
