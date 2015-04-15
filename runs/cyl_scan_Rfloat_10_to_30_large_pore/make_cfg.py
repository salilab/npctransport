#!/usr/bin/python
from IMP.npctransport import *
import sys
import math

# fetch params
# Usage: <cmd> <outfile> [kaps_R=25.0]
outfile = sys.argv[1]
kaps_R = 25.0
if(len(sys.argv) > 2):
    kaps_R = float(sys.argv[2])
print "kaps_R = %.2f" % (kaps_R)

def get_basic_config():
    config = Configuration()
    IMP.npctransport.set_default_configuration(config)
    config.statistics_fraction.lower=1.0
    #config.dump_interval=1
    config.interaction_k.lower=10
    config.interaction_range.lower=1
    # create_range(config.backbone_k, .2, 1, 10
    config.backbone_k.lower=1
    #config.time_step_factor.lower=0.3
    config.time_step_factor.lower=1
    #create_range(config.rest_length_factor, .5, 1, 10)
    config.excluded_volume_k.lower=20
    config.nonspecific_range.lower=2
    config.nonspecific_k.lower=0.01
    config.slack.lower = 8
    config.number_of_trials=1
    config.angular_D_factor.lower=0.3 #increased dynamic viscosity relative to water?
    config.statistics_interval_ns=10
    ###
    #simulation bounding volumes:
    config.box_is_on.lower=1
    return config


def add_interactions_for_fg(fg_name,
                            k_kap_lower,
                            k_kap_upper = 0, # relevant only if k_kap_steps > 1
                            k_kap_steps = 1,
                            k_kap_base = 1):
    interactionFG_KAP= IMP.npctransport.add_interaction(config,
                                     name0=fg_name,
                                     name1="kap",
                                     interaction_k=k_kap_lower,
                                     interaction_range=2)
    if(k_kap_steps > 1):
        create_range(interactionFG_KAP.interaction_k,
                     k_kap_lower, k_kap_upper,
                     steps = k_kap_steps,
                     base = k_kap_base)
    # interactionFG_CRAP= IMP.npctransport.add_interaction(config,
    #                                    name0=fg_name,
    #                                    name1="crap0",
    #                                    interaction_k=0,
    #                                    interaction_range=0)

# ********* MAIN: *********
config= get_basic_config()
#config.dump_interval_ns=0.01
#config.simulation_time_ns=5
config.dump_interval_ns=250
config.simulation_time_ns=10000
config.box_is_on.lower=1
config.box_side.lower=700
config.slab_is_on.lower=1
config.slab_thickness.lower=150
config.tunnel_radius.lower=135
#create_range(config.interaction_k, lb=.5, ub=15, steps=8, base=1)

# fg_cyto= IMP.npctransport.add_fg_type(config,
#                                  number_of_beads=12,
#                                  number=8,
#                                  radius=6,
#                                  interactions=1,
#                                  rest_length_factor = 1.5)
fg_middle= IMP.npctransport.add_fg_type(config,
                                 type_name="fg0",
                                 number_of_beads=45,
                                 number=32,
                                 radius=6,
                                 interactions=1,
                                 rest_length_factor = 1.5)
#create_range(fg_middle.number_of_beads,18,20,2)
# fg_nuclear= IMP.npctransport.add_fg_type(config,
#                                  number_of_beads=12,
#                                  number=8,
#                                  radius=6,
#                                  interactions=1,
#                                  rest_length_factor = 1.5)
#surface_area_ratio_to_R25 = math.pow( floaters_R / 25.0, 2 )
#kap_n_interactions = int( math.ceil ( 12 * surface_area_ratio_to_R25 ) )
kaps= IMP.npctransport.add_float_type(config,
                                     number=1,
                                     radius=kaps_R,
                                     interactions=12)
#create_range(kaps.radius, lb = 10, ub = 30, steps = 5, base = 1)
nonspecifics={}
rrange=range(10,30,2)
for i in rrange:
    type_name="R%d" % i
    nonspecifics[i]= IMP.npctransport.add_float_type(config,
                                                  number=50,
                                                  radius=i,
                                                  type_name=type_name,
                                                  interactions=0)
    IMP.npctransport.add_interaction(config,
                                     name0="fg0",
                                     name1=type_name,
                                     interaction_k=0,
                                     interaction_range=0)
#create_range(nonspecifics.radius, lb = 15, ub = 30, steps = 1, base = 1)
# fg with kaps / craps
add_interactions_for_fg("fg0", k_kap_lower=0.1)

#add_interactions_for_fg("fg1", kap_k)
#add_interactions_for_fg("fg2", kap_k)

# non-specific attraction
config.nonspecific_range.lower= 1.0
config.nonspecific_k.lower= 1.5

interactionFG_FG= IMP.npctransport.add_interaction(config,
                                                   name0= "fg0",
                                                   name1= "fg0",
                                                   interaction_k=0.01,
                                                   interaction_range= 2)
#create_range(interactionFG_FG.interaction_k, lb = 0.01, ub = 10, steps = 7, base = 2)
## internal FG-FG
#for i in range(3):
#    for j in range(i,3):
#        interactionFG_FG= IMP.npctransport.add_interaction(config,
#                                                           name0= "fg%d" % i,
#                                                           name1= "fg%d" % j,
#                                                           interaction_k= float(fgs_k),
#                                                           interaction_range= 2)

# dump to file
f=open(outfile, "wb")
f.write(config.SerializeToString())
print config
