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
create_range(config.time_step_factor, .001, 1000, 10)
#create_range(config.rest_length_factor, .5, 1, 10)
config.box_size.lower=100
config.slack.lower=10
config.number_of_trials.lower=1000
config.number_of_frames.lower=1000000

fg= config.fgs.add()
fg.number_of_beads.lower=10
create_range(fg.number, 2, 10, 4)
fg.radius.lower=10
fg.interactions.lower=1
fg.rest_length_factor.lower=1
fg.D_factor.lower=1

kap= config.floaters.add()
kap.number.lower=4
kap.radius.lower=20
kap.interactions.lower=10
kap.D_factor.lower=1

f=open(sys.argv[1], "wb")
f.write(config.SerializeToString())

print config
if len(sys.argv)>2:
    config.number_of_trials.lower=2
    config.number_of_frames.lower=2
    config.dump_interval=1
    f=open(sys.argv[2], "wb")
    f.write(config.SerializeToString())
