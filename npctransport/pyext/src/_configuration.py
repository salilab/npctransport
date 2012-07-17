nranges=1

def create_range(field, lb, ub=None, steps=5, base=2):
    """Create a range in the passed field"""
    field.lower=lb
    if ub:
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
    config.interaction_k.lower=1
    config.interaction_range.lower=5
    config.backbone_k.lower=1
    config.time_step_factor.lower=1
    config.box_is_on.lower=1
    config.box_side.lower=100
    config.slab_is_on.lower=0
    config.slab_thickness.lower=30
    config.tunnel_radius.lower=30
    config.slack.lower=5
    config.number_of_trials=40
    config.maximal_number_of_frames=1000000000
    config.simulation_time_nanosec=10
    config.dump_interval=10000
    config.is_dump_interval_in_ns=False
    config.nonspecific_range.lower=2
    config.nonspecific_k.lower=.3
    config.angular_D_factor.lower=5
    config.statistics_interval=1000
    config.is_statistics_interval_in_ns=False
    config.excluded_volume_k.lower=1

def add_fg_type(config, number_of_beads, number, radius,
                interactions=1, rest_length_factor=1, d_factor=1,
                interaction_k_factor=1, interaction_range_factor=1):
    fg= config.fgs.add()
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
                interaction_k_factor=1, interaction_range_factor=1):
    fg= config.floaters.add()
    fg.number.lower=number
    fg.radius.lower=radius
    fg.interactions.lower=interactions
    fg.d_factor.lower=d_factor
    fg.interaction_k_factor.lower=interaction_k_factor
    fg.interaction_range_factor.lower=interaction_range_factor
    return fg

def add_interaction(config, name0, name1,
                    interaction_k=None, interaction_range=None, is_on=1):
    i= config.interactions.add()
    i.type0= name0
    i.type1=name1
    if interaction_k:
      i.interaction_k.lower=interaction_k
    if interaction_range:
      i.interaction_range.lower=interaction_range
    i.is_on.lower=is_on
    return i


def set_single_configuration(config):
    """Change the passed configuration to be a single run configuration"""
    config.number_of_trials=1

def set_quick_configuration(config):
    """Change the passed configuration to be quick"""
    config.maximal_number_of_frames=100000

def write(config):
    import sys
    import optparse
    parser = optparse.OptionParser(usage="usage: %prog [options] output.pb")
    parser.add_option("-s", "--single", dest="single",
                      help="Where to put a protobuf to do a single run")
    parser.add_option("-q", "--quick", dest="quick",
                      help="Where to put the protobuf for a single quick run")
    (options, args) = parser.parse_args()
    if len(args) != 1:
        parser.print_help()
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
    print nranges, "work units"
