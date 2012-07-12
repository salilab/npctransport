#!/usr/bin/python
from IMP.npctransport import *
import sys
import barak_basic_configuration

config= barak_basic_configuration.get_basic_config()
config.dump_interval=50000
config.number_of_frames=200000000

fg= IMP.npctransport.add_fg_type(config,
                                 number_of_beads=12,
                                 number=12,
                                 radius=8,
                                 interactions=1)
kap= IMP.npctransport.add_float_type(config,
                                     number=10,
                                     radius=20,
                                     interactions=12)
#nonspecifics= IMP.npctransport.add_float_type(config,
#                                              number=10,
#                                              radius=20,
#                                              interactions=0)

interactionFG_KAP= IMP.npctransport.add_interaction(config,
                                                    name0="fg0",
                                                    name1="kap",
                                                    interaction_k=30,
                                                    interaction_range=5)
#interactionFG_CRAP= IMP.npctransport.add_interaction(config,
#                                                    name0="fg0",
#                                                    name1="crap0",
#                                                    interaction_k=0,
#                                                    interaction_range=0)
interactionFG_FG= IMP.npctransport.add_interaction(config,
                                                   name0="fg0",
                                                   name1="fg0",
                                                   interaction_k=30,
                                                   interaction_range=1)
#interactionCRAP_KAP= IMP.npctransport.add_interaction(config,
#                                                   name0="crap0",
#                                                   name1="kap",
#                                                   interaction_k=0,
#                                                   interaction_range=0)

create_range(config.nonspecific_k, 0.01, 1, 5)
create_range(config.nonspecific_range, 0.1, 10, 5)
create_range(interactionFG_KAP.interaction_k, 0.1, 500, 50)
create_range(interactionFG_KAP.interaction_range, 0.1, 10, 10)
create_range(interactionFG_FG.interaction_k, 0.1, 500, 30)
create_range(interactionFG_FG.interaction_range, 0.1, 10, 10)


f=open(sys.argv[1], "wb")
f.write(config.SerializeToString())

print config
if len(sys.argv)>2:
    config.number_of_trials.lower=2
    config.number_of_frames.lower=2
    config.dump_interval=1
    f=open(sys.argv[2], "wb")
    f.write(config.SerializeToString())
