#!/usr/bin/python
from IMP.npctransport import *
import sys
def create_range(field, lb, ub=None, steps=None, base=None):
    field.lower=lb
    if ub:
        field.upper=ub
        field.steps=steps
        if base:
            field.base=base
        else:
            try:
                field.base=2
            except:
                pass

config= Configuration()
config.dump_interval=400
create_range(config.interaction_spring_constant, .2, 20, 10)
create_range(config.backbone_spring_constant, .2, 20, 10)
config.time_step_factor.lower=1
create_range(config.nonspecific_range, 0, 5, 3, base=1)
#create_range(config.rest_length_factor, .5, 1, 10)
config.box_size.lower=1000
config.slack.lower=10
config.number_of_trials.lower=1000
config.number_of_frames.lower=1000

fg= config.fgs.add()
create_range(fg.number_of_beads,2, 20,10)
fg.number.lower=1
create_range(fg.radius, 5,50, 5, base=1)
fg.interactions.lower=1
create_range(fg.rest_length_factor, .5, 1, 3)
fg.D_factor.lower=1

interaction= config.interactions.add()
interaction.type0="fg0"
interaction.type1="fg0"
create_range(interaction.on_or_off, 0, 1, 2)


f=open(sys.argv[1], "wb")
f.write(config.SerializeToString())

print config
if len(sys.argv)>2:
    config.number_of_trials.lower=1
    config.number_of_frames.lower=1
    config.dump_interval=1
    f=open(sys.argv[2], "wb")
    f.write(config.SerializeToString())
