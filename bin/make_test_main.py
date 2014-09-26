#!/usr/bin/env python
from IMP.npctransport import *
import sys
"""Create the protobuf for the test test_main"""

nonspecifics_on=False

config= Configuration()
set_default_configuration(config)
config.dump_interval_ns=0.5 # in ns
config.statistics_interval_ns=0.01 # in ns
config.maximal_number_of_frames=100
config.simulation_time_ns=.00001
config.box_side.lower = 0.3 * config.box_side.lower
config.time_step_factor.lower=3

fg= IMP.npctransport.add_fg_type(config,
                                 number_of_beads=2,
                                 number=1,
                                 radius=8,
                                 interactions=1)
kap= IMP.npctransport.add_float_type(config,
                                     number=1,
                                     radius=20,
                                     interactions=12)
if(nonspecifics_on):
    nonspecifics= IMP.npctransport.add_float_type(config,
                                                  number=2,
                                                  radius=20,
                                                  interactions=0)
    interactionFG_CRAP= IMP.npctransport.add_interaction(config,
                                                         name0="fg0",
                                                         name1="crap0",
                                                         interaction_k=0,
                                                         interaction_range=0)

interactionFG_KAP= IMP.npctransport.add_interaction(config,
                                                    name0="fg0",
                                                    name1="kap",
                                                    interaction_k=30,
                                                    interaction_range=2)
interactionFG_FG= IMP.npctransport.add_interaction(config,
                                                   name0="fg0",
                                                   name1="fg0",
                                                   interaction_k=30,
                                                   interaction_range=2)

#interactionFG_KAP.interaction_k.lower = 2000
create_range(interactionFG_KAP.interaction_k, 0.001, 500, steps=30,base=1.3)
#create_range(interactionFG_KAP.interaction_range, 0.1, 6, steps=5)
create_range(config.nonspecific_k, 0.0001, 1, steps=10)
#create_range(config.nonspecific_range, 0.1, 5, steps=5)
#create_range(interactionFG_FG.interaction_k, 0.001, 500, steps=20)
#create_range(interactionFG_FG.interaction_range, 0.1, 6, steps=5)

f=open(sys.argv[1], "wb")
f.write(config.SerializeToString())

print config
if len(sys.argv)>2:
    config.number_of_trials.lower=2
    config.number_of_frames.lower=2
    config.dump_interval=1
    f=open(sys.argv[2], "wb")
    f.write(config.SerializeToString())
