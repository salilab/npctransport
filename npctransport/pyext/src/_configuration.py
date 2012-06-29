def create_range(field, lb, ub=None, steps=None, base=None):
    """Create a range in the passed field"""
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

def set_default_configuration(config):
    """Set the defaults for the configuration"""
    config.interaction_k.lower=1
    config.interaction_range.lower=5
    config.backbone_k.lower=1
    config.time_step_factor.lower=1
    config.box_on_or_off.lower=1
    config.box_side.lower=100
    config.slab_on_or_off.lower=0
    config.slab_thickness.lower=30
    config.tunnel_radius.lower=30
    config.slack.lower=5
    config.number_of_trials=40
    config.number_of_frames=1000000
    config.dump_interval=10000
    config.nonspecific_range.lower=2
    config.nonspecific_k.lower=.3
    config.angular_D_factor.lower=5
    config.statistics_interval=1000
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

def add_interaction(config, name0, name1, is_on=1):
    i= config.interactions.add()
    i.type0= name0
    i.type1=name1
    i.is_on.lower=is_on
    return i


def set_single_configuration(config):
    """Change the passed configuration to be a single run configuration"""
    config.number_of_trials.lower=1

def set_quick_configuration(config):
    """Change the passed configuration to be quick"""
    config.number_of_frames.lower=100000

def write(config):
    import sys
    f=open(sys.argv[1], "wb")
    f.write(config.SerializeToString())

    if len(sys.argv)>2:
        set_single_configuration(config)
        f=open(sys.argv[2], "wb")
        f.write(config.SerializeToString())
    if len(sys.argv)>3:
        set_quick_configuration(config)
        #config.dump_interval=1
        f=open(sys.argv[3], "wb")
        f.write(config.SerializeToString())
