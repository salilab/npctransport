#!/usr/bin/python
from IMP.npctransport import *
import sys
import barak_basic_configuration

nonspecifics_on=None
if len(sys.argv)>2:
    nonspecifics_on = (sys.argv[2] != '0')
    print "Nonspecifics: ", nonspecifics_on
else:
    print "Usage: %s <output-file> <is_nonspecifics_on> [output-file-quick]" \
        % sys.argv[0]

config= barak_basic_configuration.get_basic_config()
config.simulation_time_ns=5000
config.dump_interval_ns=500
config.statistics_interval_ns=0.1
config.maximal_number_of_frames=400000000
config.time_step_factor.lower=2.5

fg= IMP.npctransport.add_fg_type(config,
                                 number_of_beads=12,
                                 number=12,
                                 radius=8,
                                 interactions=1)
kap= IMP.npctransport.add_float_type(config,
                                     number=10,
                                     radius=20,
                                     interactions=12)
if(nonspecifics_on):
    nonspecifics= IMP.npctransport.add_float_type(config,
                                                  number=10,
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
                                                    interaction_range=5)
interactionFG_FG= IMP.npctransport.add_interaction(config,
                                                   name0="fg0",
                                                   name1="fg0",
                                                   interaction_k=30,
                                                   interaction_range=1)

create_range(interactionFG_KAP.interaction_k, 0.0001, 10, steps=30)
create_range(interactionFG_KAP.interaction_range, 0.1, 6, steps=5)
create_range(config.nonspecific_k, 0.0001, 10, steps=10)
create_range(config.nonspecific_range, 0.1, 6, steps=5)
create_range(interactionFG_FG.interaction_k, 0.0001, 30, steps=10)
create_range(interactionFG_FG.interaction_range, 0.1, 6, steps=5)
"""
create_range(interactionFG_KAP.interaction_k, 0.0001, 5, steps=20)
create_range(interactionFG_KAP.interaction_range, 0.5, 5, steps=3)
create_range(config.nonspecific_k, 0.0001, 100, steps=20)
create_range(config.nonspecific_range, 0.1, 6, steps=5)
create_range(interactionFG_FG.interaction_k, 0.0001, 30, steps=10)
create_range(interactionFG_FG.interaction_range, 0.5, 5, steps=3)
"""


f=open(sys.argv[1], "wb")
f.write(config.SerializeToString())

print config
if len(sys.argv)>3:
    config.number_of_trials.lower=2
    config.number_of_frames.lower=2
    config.dump_interval=1
    f=open(sys.argv[2], "wb")
    f.write(config.SerializeToString())
