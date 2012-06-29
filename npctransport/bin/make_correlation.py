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
config.dump_interval=1000
config.statistics_interval=100
create_range(config.interaction_k, .2, 20, 5)
create_range(config.backbone_k, .2, 20, 5)
config.time_step_factor.lower=1
create_range(config.nonspecific_range, 0, 5, 3, base=1)
create_range(config.nonspecific_k, .2, 20, 3, 5)
create_range(config.angular_D_factor, 10, 50, 3)
#create_range(config.rest_length_factor, .5, 1, 10)
config.box_size.lower=150
config.slack.lower=10
config.number_of_trials.lower=10
config.number_of_frames.lower=1000000

fg= config.fgs.add()
fg.number_of_beads.lower=10
fg.number.lower=1
fg.radius.lower=10
fg.interactions.lower=1
create_range(fg.rest_length_factor, .5, .8, 3)
fg.D_factor.lower=1

kap= config.floaters.add()
create_range(kap.number,0, 10,3)
kap.radius.lower=10
create_range(kap.interactions, 1, 10, 3)

interaction= config.interactions.add()
interaction.type0="fg0"
interaction.type1="fg0"
create_range(interaction.on_or_off, 0, 1, 2)

interaction= config.interactions.add()
interaction.type0="fg0"
interaction.type1="kap"
create_range(interaction.on_or_off, 0, 1, 2)


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
