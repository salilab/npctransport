#!/usr/bin/env python

import sys
import IMP.npctransport


config= IMP.npctransport.Configuration()
IMP.npctransport.set_default_configuration(config)
IMP.npctransport.create_range(config.interaction_k, .2, 20, 5)
IMP.npctransport.create_range(config.backbone_k, .2, 20, 5)
IMP.npctransport.create_range(config.nonspecific_range, 0, 5, 3, base=1)
IMP.npctransport.create_range(config.nonspecific_k, .2, 20, 3)
config.box_side.lower=400
config.number_of_trials=1
config.number_of_frames=10000
config.dump_interval=100
config.slab_thickness.lower=30;
config.tunnel_radius.lower=20;
config.slab_is_on.lower=1;

#IMP.npctransport.create_range(config.angular_D_factor, 10, 50, 3)

fg= IMP.npctransport.add_fg_type(config,
                                 number_of_beads=10,
                                 number=40,
                                 radius=10,
                                 interactions=1)
IMP.npctransport.create_range(fg.rest_length_factor, .5, .8, 3)

kap= IMP.npctransport.add_float_type(config,
                                  number=40,
                                  radius=10)
kap.interactions.lower=10
crap0= IMP.npctransport.add_float_type(config,
                                      number=30,
                                      radius=10)
crap0.interactions.lower=0
crap1= IMP.npctransport.add_float_type(config,
                                      number=30,
                                      radius=20)
crap1.interactions.lower=1

interaction= IMP.npctransport.add_interaction(config, "fg0", "fg0")
IMP.npctransport.create_range(interaction.is_on,0,1,2)

interaction= IMP.npctransport.add_interaction(config, "fg0", "kap")

config.statistics_interval=1000000000

IMP.npctransport.write(config)
