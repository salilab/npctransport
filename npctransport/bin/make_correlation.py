#!/usr/bin/python

import sys


config= Configuration()
IMP.npctransport.set_default_configuration(config)
IMP.npctransport.create_range(config.interaction_k, .2, 20, 5)
IMP.npctransport.create_range(config.backbone_k, .2, 20, 5)
IMP.npctransport.create_range(config.nonspecific_range, 0, 5, 3, base=1)
IMP.npctransport.create_range(config.nonspecific_k, .2, 20, 3, 5)
IMP.npctransport.create_range(config.angular_D_factor, 10, 50, 3)

fg= IMP.npctransport.add_fg_type(config,
                                 number_of_beads=10,
                                 number=1,
                                 radius=10,
                                 interactions=1)
IMP.npctransport.create_range(fg.rest_length_factor, .5, .8. 3)

kap= IMP.npctransport.add_floater(config,
                                  number_of_beads=10,
                                  number=1,
                                  radius=10)
IMP.npctransport.create_range(kap.number,0, 10,3)
IMP.npctransport.create_range(kap.interactions, 1, 10, 3)

interaction= IMP.npctransport.add_interaction(config, "fg0", "fg0")
IMP.npctransport.create_range(interaction.is_on,0,1,2)

interaction= IMP.npctransport.add_interaction(config, "fg0", "kap")

f=open(sys.argv[1], "wb")
f.write(config.SerializeToString())

print config
if len(sys.argv)>2:
    config.number_of_trials.lower=1
    f=open(sys.argv[2], "wb")
    f.write(config.SerializeToString())
if len(sys.argv)>3:
    config.number_of_frames.lower=100000
    #config.dump_interval=1
    f=open(sys.argv[3], "wb")
    f.write(config.SerializeToString())
