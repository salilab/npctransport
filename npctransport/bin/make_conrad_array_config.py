#!/usr/bin/python
from IMP.npctransport import *
import sys
def create_range(field, lb):
    field.lower=lb

config= Configuration()
config.dump_interval=400
create_range(config.interaction_spring_constant, 1)
create_range(config.backbone_spring_constant, 1)
config.time_step_factor.lower=10
create_range(config.nonspecific_range, 1)
#create_range(config.rest_length_factor, .5, 1, 10)
config.box_size.lower=200
config.slack.lower=10
config.number_of_trials.lower=1
config.number_of_frames.lower=10000

fg= config.fgs.add()
create_range(fg.number_of_beads,8)
create_range(fg.number, 10)
create_range(fg.radius, 5)
fg.interactions.lower=1
create_range(fg.rest_length_factor, .5)
config.angular_D_factor.lower=10

kap= config.floaters.add()
create_range(kap.number,8)
create_range(kap.radius, 5)
kap.interactions.lower=10

interaction= config.interactions.add()
interaction.type0="fg0"
interaction.type1="fg0"
create_range(interaction.on_or_off, 1)

interaction= config.interactions.add()
interaction.type0="fg0"
interaction.type1="kap"
create_range(interaction.on_or_off, 1)

config.dump_interval=10

f=open(sys.argv[1], "wb")
f.write(config.SerializeToString())

print config
