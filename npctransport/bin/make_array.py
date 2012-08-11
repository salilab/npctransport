#!/usr/bin/python

import sys
import IMP.npctransport

# 121500

config= IMP.npctransport.Configuration()
IMP.npctransport.set_default_configuration(config)
IMP.npctransport.create_range(config.interaction_k, .2, 20, 5)
IMP.npctransport.create_range(config.interaction_range, 1, 5, 3)
IMP.npctransport.create_range(config.backbone_k, .2, 20, 5)
IMP.npctransport.create_range(config.nonspecific_range, 0, 5, 3, base=1)
IMP.npctransport.create_range(config.nonspecific_k, .2, 20, 3)

config.box_side.lower=150
config.simulation_time_ns=3000
config.number_of_trials=4
config.dump_interval_ns=10

#IMP.npctransport.create_range(config.angular_D_factor, 10, 50, 3)

fg= IMP.npctransport.add_fg_type(config,
                                 number_of_beads=10,
                                 number=10,
                                 radius=10,
                                 interactions=1)
kap= IMP.npctransport.add_float_type(config,
                                  number=4,
                                  radius=10)
IMP.npctransport.create_range(kap.interactions, 2, 10, 5)
IMP.npctransport.create_range(kap.number, 2, 10, 3)
crap0= IMP.npctransport.add_float_type(config,
                                       number=30,
                                       radius=10)
IMP.npctransport.create_range(crap0.number, 0, 100, 6)
crap0.interactions.lower=1

interaction= IMP.npctransport.add_interaction(config, "fg0", "fg0")
IMP.npctransport.create_range(interaction.is_on,0,1,2)

interaction= IMP.npctransport.add_interaction(config, "fg0", "kap")
interaction= IMP.npctransport.add_interaction(config, "fg0", "crap0")

config.statistics_interval_ns=1

IMP.npctransport.write(config)
