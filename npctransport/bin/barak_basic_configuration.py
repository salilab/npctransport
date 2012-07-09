#!/usr/bin/python
from IMP.npctransport import *
import sys

def get_basic_config():
    config = Configuration()
    IMP.npctransport.set_default_configuration(config)
    config.dump_interval=50000
    #config.dump_interval=1
    config.interaction_k.lower=10
    config.interaction_range.lower=1
    # create_range(config.backbone_k, .2, 1, 10)
    config.backbone_k.lower=.02
    #config.time_step_factor.lower=0.3
    config.time_step_factor.lower=3
    #create_range(config.rest_length_factor, .5, 1, 10)
    config.excluded_volume_k.lower=.02
    config.nonspecific_range.lower=2
    config.nonspecific_k.lower=0.01
    create_range(config.slack, 10, 11, 2)
    config.number_of_trials=1
    config.number_of_frames=200000000
    #config.number_of_frames.lower=500
    config.angular_D_factor.lower=200
    config.statistics_interval=1000
    ###
    #simulation bounding volumes:
    config.box_is_on.lower=1
    config.box_side.lower=350
    config.slab_is_on.lower=0
    config.slab_thickness.lower=225
    config.tunnel_radius.lower=100
    return config
