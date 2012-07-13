#!/usr/bin/python
from IMP.npctransport import *
import sys
import barak_basic_configuration

def add_interactions_for_fg(fg_name, k_kap):
    interactionFG_KAP= IMP.npctransport.add_interaction(config,
                                                        name0=fg_name,
                                                        name1="kap",
                                                        interaction_k=k_kap,
                                                        interaction_range=3)
    interactionFG_CRAP= IMP.npctransport.add_interaction(config,
                                                         name0=fg_name,
                                                         name1="crap0",
                                                         interaction_k=0,
                                                         interaction_range=0)
    interactionFG_FG= IMP.npctransport.add_interaction(config,
                                                       name0=fg_name,
                                                       name1=fg_name,
                                                       interaction_k=30,
                                                       interaction_range=1)

# ********* MAIN: *********
config= barak_basic_configuration.get_basic_config()
config.dump_interval=25000
config.simulation_time_nanosec=500
config.box_is_on.lower=1
config.box_side.lower=350
config.slab_is_on.lower=1
config.slab_thickness.lower=225
config.tunnel_radius.lower=100

fg_cyto= IMP.npctransport.add_fg_type(config,
                                 number_of_beads=8,
                                 number=6,
                                 radius=8,
                                 interactions=1)
fg_middle= IMP.npctransport.add_fg_type(config,
                                 number_of_beads=8,
                                 number=12,
                                 radius=8,
                                 interactions=1)
fg_nuclear= IMP.npctransport.add_fg_type(config,
                                 number_of_beads=12,
                                 number=6,
                                 radius=8,
                                 interactions=1)
kap= IMP.npctransport.add_float_type(config,
                                     number=10,
                                     radius=20,
                                     interactions=12)
nonspecifics= IMP.npctransport.add_float_type(config,
                                              number=10,
                                              radius=20,
                                              interactions=0)

add_interactions_for_fg("fg0",90)
add_interactions_for_fg("fg1",10)
add_interactions_for_fg("fg2",90)


f=open(sys.argv[1], "wb")
f.write(config.SerializeToString())

print config
if len(sys.argv)>2:
    config.number_of_trials.lower=2
    config.number_of_frames.lower=2
    config.dump_interval=1
    f=open(sys.argv[2], "wb")
    f.write(config.SerializeToString())
