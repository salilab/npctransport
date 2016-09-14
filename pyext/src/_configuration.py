from __future__ import print_function

nranges=1

def create_range(field, lb, ub=None, steps=5, base=2):
    """
    Create a range in the passed field
    using a logarithmic-evenely distributed steps
    (or even-sized steps for base=1)

    lb - lower bound value
    ub - upper bound value, if None then only lower bound is used
    steps - number of steps from lb to ub, if 1 then only lower bound is used
    base - lograithmic base for spacing, use 1 for standard even-sized steps
    """
    field.lower=lb
    if ub and steps != 1:
        field.upper=ub
        field.steps=steps
        global nranges
        nranges*=steps
        try:
            field.base=base
        except:
            # int range
            pass


def set_default_configuration(config):
    """Set the defaults for the configuration"""
    default_version=config.output_npctransport_version
    config.output_npctransport_version=default_version
    config.interaction_k.lower=1
    config.interaction_range.lower=5
    config.backbone_k.lower=5
    config.time_step_factor.lower=1
    config.box_is_on.lower=1
    config.box_side.lower=100
    config.slab_is_on.lower=0
    config.slab_thickness.lower=30
    config.tunnel_radius.lower=30
    config.slack.lower=15
    config.number_of_trials=40
    config.maximal_number_of_frames=1000000000
    config.simulation_time_ns=10
    config.dump_interval_ns=50
    config.nonspecific_range.lower=8
    config.nonspecific_k.lower=.075
    config.angular_D_factor.lower=5
    config.statistics_interval_ns=.001
    config.excluded_volume_k.lower=20
    config.statistics_fraction.lower=.5;
    config.fg_anchor_inflate_factor=1.0
    config.are_floaters_on_one_slab_side=0; # false
    config.is_exclude_floaters_from_slab_initially=1; # true

def add_fg_type(config, type_name, number_of_beads, number, radius,
                interactions=1, rest_length_factor=1, d_factor=1,
                interaction_k_factor=1, interaction_range_factor=1):
    fg= config.fgs.add()
    fg.type = type_name
    fg.number_of_beads.lower=number_of_beads
    fg.number.lower=number
    fg.radius.lower=radius
    fg.interactions.lower=interactions
    fg.rest_length_factor.lower=rest_length_factor
    fg.d_factor.lower=d_factor
    fg.interaction_k_factor.lower=interaction_k_factor
    fg.interaction_range_factor.lower=interaction_range_factor
    return fg

def add_float_type(config, number, radius,
                   interactions=1,  d_factor=1,
                interaction_k_factor=1, interaction_range_factor=1,
                   type_name=None):
    f= config.floaters.add()
    f.number.lower=number
    f.radius.lower=radius
    f.interactions.lower=interactions
    f.d_factor.lower=d_factor
    f.interaction_k_factor.lower=interaction_k_factor
    f.interaction_range_factor.lower=interaction_range_factor
    if type_name != None:
        f.type = type_name
    return f

def add_obstacle_type(config, type_name, R, xyzs=[],
                      is_static=1,
                      interactions=0,  d_factor=1,
                      interaction_k_factor=1, interaction_range_factor=1):
    """
    add an obstacle type name type_name, with radius R
    Note that the .xyzs vector field must be filled in to get actual instances
    and is by default empty

    params:
    xyzs - an array of x,y,z coordinates triplets
    """
    o= config.obstacles.add()
    o.type = type_name
    o.radius.lower= R
    for xyz in xyzs:
        assert(len(xyz) == 3)
        new_xyz=o.xyzs.add()
        new_xyz.x = xyz[0]
        new_xyz.y = xyz[1]
        new_xyz.z = xyz[2]
    o.interactions.lower= interactions
    o.is_static= is_static
    o.d_factor.lower= d_factor
    o.interaction_k_factor.lower= interaction_k_factor
    o.interaction_range_factor.lower= interaction_range_factor
    return o

def add_interaction(config, name0, name1,
                    interaction_k=None, interaction_range=None, is_on=1,
                    range_sigma0_deg=None, range_sigma1_deg=None):
    i= config.interactions.add()
    i.type0= name0
    i.type1=name1
    if interaction_k:
      i.interaction_k.lower=interaction_k
    if interaction_range:
      i.interaction_range.lower=interaction_range
    i.is_on.lower=is_on
    if not (range_sigma0_deg is None and range_sigma1_deg is None):
        i.range_sigma0_deg.lower = range_sigma0_deg
        i.range_sigma1_deg.lower = range_sigma1_deg
    return i


def set_single_configuration(config):
    """Change the passed configuration to be a single run configuration"""
    config.number_of_trials=1

def set_quick_configuration(config):
    """Change the passed configuration to be quick"""
    config.maximal_number_of_frames=100000

import sys
import optparse
make_parser = optparse.OptionParser(usage="usage: %prog [options] output.pb")
make_parser.add_option("-s", "--single", dest="single",
                  help="Where to put a protobuf to do a single run")
make_parser.add_option("-q", "--quick", dest="quick",
                  help="Where to put the protobuf for a single quick run")

def write(config):
    (options, args) = make_parser.parse_args()
    if len(args) != 1:
        make_parser.print_help()
        exit(1)
    f=open(args[0], "wb")
    f.write(config.SerializeToString())

    if options.single:
        set_single_configuration(config)
        f=open(options.single, "wb")
        f.write(config.SerializeToString())
    if options.quick:
        set_quick_configuration(config)
        #config.dump_interval=1
        f=open(options.quick, "wb")
        f.write(config.SerializeToString())
    print(nranges, "work units")
