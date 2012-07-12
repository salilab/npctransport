#!/usr/bin/python
from IMP.npctransport import *
import sys
import barak_basic_configuration

config= barak_basic_configuration.get_basic_config()
config.dump_interval=25000
config.number_of_frames=50000000
config.box_is_on.lower=1
config.box_side.lower=350
config.slab_is_on.lower=1
config.slab_thickness.lower=225
config.tunnel_radius.lower=100

fg= IMP.npctransport.add_fg_type(config,
                                 number_of_beads=8,
                                 number=24,
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

interactionFG_KAP= IMP.npctransport.add_interaction(config,
                                                    name0="fg0",
                                                    name1="kap",
                                                    interaction_k=90,
                                                    interaction_range=3)
interactionFG_CRAP= IMP.npctransport.add_interaction(config,
                                                     name0="fg0",
                                                     name1="crap0",
                                                     interaction_k=0,
                                                     interaction_range=0)
interactionFG_FG= IMP.npctransport.add_interaction(config,
                                                   name0="fg0",
                                                   name1="fg0",
                                                   interaction_k=30,
                                                   interaction_range=1)
interactionCRAP_KAP= IMP.npctransport.add_interaction(config,
                                                      name0="crap0",
                                                      name1="kap",
                                                      interaction_k=0,
                                                      interaction_range=0)


f=open(sys.argv[1], "wb")
f.write(config.SerializeToString())

print config
if len(sys.argv)>2:
    config.number_of_trials.lower=2
    config.number_of_frames.lower=2
    config.dump_interval=1
    f=open(sys.argv[2], "wb")
    f.write(config.SerializeToString())
