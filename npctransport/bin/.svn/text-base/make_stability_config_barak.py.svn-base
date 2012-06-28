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
config.interaction_spring_constant.lower=10
# create_range(config.backbone_spring_constant, .2, 1, 10)
config.backbone_spring_constant.lower=.02
config.time_step_factor.lower=0.3
#create_range(config.rest_length_factor, .5, 1, 10)
config.nonspecific_range.lower=1
config.slack.lower=10
config.number_of_trials.lower=1
config.number_of_frames.lower=1000000
config.angular_D_factor.lower=200
config.statistics_interval=100

#simulation bounding volumes:
config.box.on_or_off.lower=0
config.box.size.lower=300
config.slab.on_or_off.lower=1
config.slab.height.lower=600
config.slab.radius.lower=100
config.slab.width.lower=350


fg= config.fgs.add()
fg.number_of_beads.lower=10
fg.number.lower=25
fg.radius.lower=8
fg.interactions.lower=1
fg.rest_length_factor.lower=1
fg.D_factor.lower=1

kap=config.floaters.add()
kap.number.lower=10
kap.radius.lower=20
kap.interactions.lower=10
kap.D_factor.lower=1

nonspecifics= config.floaters.add()
nonspecifics.number.lower=10
nonspecifics.radius.lower=20
#create_range( nonspecifics.radius, 10, 30, 10 )
nonspecifics.interactions.lower=0
nonspecifics.D_factor.lower=1

interactionFG_KAP = config.interactions.add()
interactionFG_KAP.type0="fg0"
interactionFG_KAP.type1="kap"
interactionFG_KAP.on_or_off.lower=1

interactionFG_FG = config.interactions.add()
interactionFG_FG.type0="fg0"
interactionFG_FG.type1="fg0"
interactionFG_FG.on_or_off.lower=1

f=open(sys.argv[1], "wb")
f.write(config.SerializeToString())

print config
if len(sys.argv)>2:
    config.number_of_trials.lower=2
    config.number_of_frames.lower=2
    config.dump_interval=1
    f=open(sys.argv[2], "wb")
    f.write(config.SerializeToString())
