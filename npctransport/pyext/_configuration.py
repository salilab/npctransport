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
    config.backbone_k=1
    config.time_step_factor=1
    config.box_on_or_off.lower=1
    config.box_size=100
    config.slab_on_or_off.lower=0
    config.slab_thickness=30
    config.tunnel_radius=30
    config.slack.lower=5
    config.number_of_trials.lower=40
    config.number_of_frames.lower=1000000
    config.dump_interval=10000
    config.nonspecific_range.lower=2
    config.nonspecific_k.lower=.3
    config.angular_D_factor.lower=5
    config.statistics_interval=1000
    config.excluded_volume_k.lower=1
    config.slab.on_or_off=0


def set_single_configuration(config):
    """Change the passed configuration to be a single run configuration"""
    config.number_of_trials.lower=1

def set_quick_configuration(config):
    """Change the passed configuration to be quick"""
    config.number_of_frames.lower=100000
