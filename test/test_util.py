from __future__ import print_function
from IMP.npctransport import *
import sys
import math
import numpy as np

def write_config_file(outfile, config):
    with open(outfile, "wb") as f:
        f.write(config.SerializeToString())

def get_basic_config():
    config = Configuration()
    IMP.npctransport.set_default_configuration(config)
    config.statistics_fraction.lower=0.9
    #config.dump_interval=1
    config.interaction_k.lower=10
    config.interaction_range.lower=1
    # create_range(config.backbone_k, .2, 1, 10
    config.backbone_k.lower=1
    #config.time_step_factor.lower=0.3
    config.time_step_factor.lower=3
    #create_range(config.rest_length_factor, .5, 1, 10)
    config.excluded_volume_k.lower=10
    config.nonspecific_range.lower=5
    config.nonspecific_k.lower=0.01
    config.slack.lower = 8
    config.number_of_trials=1
    config.dump_interval_ns=0.1
    config.simulation_time_ns=500
    config.angular_D_factor.lower=0.3 #increased dynamic viscosity relative to water?
    config.statistics_interval_ns=0.01
    ###
    #simulation bounding volumes:
    config.box_is_on.lower=1
    config.box_side.lower=400
    config.slab_is_on.lower=0
    config.slab_thickness.lower=150
    config.tunnel_radius.lower=75
    return config

def make_simple_cfg(outfile=None, is_slab_on = True, n_particles_factor = 1, is_obstacles=False,
                    fg_name = "my_fg", kap_name = "kap0", inert_name="inert0", obstacle_name="my_obstacle0"):
    """
    Make a simple configuration, with or without tunnen-in-a-slab, with one or more
    FG chains that interact with one or more kaps and inert particles, and with one
    or more static obstables, all enclosed in a bounding box and optionally a porus slab.

    Params:
    -------
    outfile - output file for configuration, if not None
    is_slab_on - suggestive
    n_particles_factor - factor by which the number of particles in the system
                         is multiplied, to generate a more complex configurations
    is_obstacles - whether to include static particles
    fg_name - name of fg type
    kap_name - name of NTRs type (facilitated diffusion)
    inert_name - name of passively diffusing molecules type
    obstacle_name - name of static obstacles (e.g. NPC scaffold in production runs)
    """
    def add_interactions_for_fg(fg_name,
                                k_kap_lower,
                                k_kap_upper = 0, # relevant only if k_kap_steps > 1
                                k_kap_steps = 1):
        """
        add interaction between fg_name and kaps / craps, with specified
        k_kap range
        """
        interactionFG_KAP= IMP.npctransport.add_interaction(config,
                                                            name0=fg_name,
                                                            name1=kap_name,
                                                            interaction_k=k_kap_lower,
                                                            interaction_range=9)
        if(k_kap_steps > 1):
            create_range(interactionFG_KAP.interaction_k,
                         k_kap_lower, k_kap_upper,
                         steps = k_kap_steps)
        interactionFG_CRAP= IMP.npctransport.add_interaction(config,
                                                             name0=fg_name,
                                                             name1=inert_name,
                                                             interaction_k=0,
                                                             interaction_range=0)

    # ********* MAIN make_simple_cfg: *********
    config= get_basic_config()
    config.dump_interval_ns=1
    config.simulation_time_ns=0.01
    config.box_is_on.lower=1
    config.box_side.lower=200
    config.slab_is_on.lower= is_slab_on
    config.slab_thickness.lower=150
    config.tunnel_radius.lower=90
    fg= IMP.npctransport.add_fg_type(config,
                                 type_name=fg_name,
                                 number_of_beads=int(math.floor(2 * n_particles_factor)),
                                 number=int(math.floor(1 * n_particles_factor)),
                                 radius=6,
                                 interactions=1,
                                 rest_length_factor = 1.5)
    kap= IMP.npctransport.add_float_type(config,
                                         type_name=kap_name,
                                         number=int(math.ceil(1 * n_particles_factor)),
                                         radius=20,
                                         interactions=6)
    nonspecifics= IMP.npctransport.add_float_type(config,
                                                  type_name=inert_name,
                                                  number=int(math.ceil(1 * n_particles_factor)),
                                                  radius=20,
                                                  interactions=0)
    if is_obstacles:
        obstacle_xyzs = [[10.0,30.0,0.0], [-10.0,30.0,0.0]]
        obstacle= IMP.npctransport.add_obstacle_type(config,
                                                     type_name=obstacle_name,
                                                     R=20,
                                                     xyzs = obstacle_xyzs)
    ###########
    # fg with kaps / craps
    add_interactions_for_fg(fg_name,
                            k_kap_lower=0.1)
     #############

    # non-specific attraction
    config.nonspecific_range.lower= 5.0
    config.nonspecific_k.lower= 0.01

    # internal FG-FG
    interactionFG_FG= IMP.npctransport.add_interaction(config,
                                                       name0= fg_name,
                                                       name1= fg_name,
                                                       interaction_k= 0.1,
                                                       interaction_range= 6)
    ##############

    # dump to file
    if outfile is not None:
        write_config_file(outfile, config)
    return config


def create_diffusing_rb_particle(m, radius):
    '''
    create a diffusing rigid-body particle of specified radius
    and mass 1.0, with hierarchy traits
    '''
    p= IMP.Particle(m)
    d= IMP.core.XYZR.setup_particle(p)
    d.set_radius(radius)
    IMP.atom.Hierarchy.setup_particle(p)
    IMP.atom.Mass.setup_particle(p, 1.0)
    IMP.core.RigidBody.setup_particle(p, IMP.algebra.ReferenceFrame3D())
    d.set_coordinates_are_optimized(True)
    IMP.atom.RigidBodyDiffusion.setup_particle(p)
    return p

def create_rb(m, radius):
    '''
    create a rigid-body particle of specified radius
    returns rigid body decorator
    '''
    p= IMP.Particle(m)
    xyzr= IMP.core.XYZR.setup_particle(p)
    xyzr.set_radius(radius)
    rb= IMP.core.RigidBody.setup_particle(p, IMP.algebra.ReferenceFrame3D())
    rb.set_coordinates_are_optimized(True)
    return rb

def test_protobuf_installed(test_class):
    '''
    Test that protobuf is installed on python. test_class is a unittest class. For typical tests, use 'self'
    '''
    protobuf_installed=False
    try:
        import google.protobuf
        protobuf_installed=True
    except ImportError:
        msg='ERROR: npctransport python module requires the python protobuf package.' \
            ' One common way is using pip - "pip install protobuf" on a local' \
            ' python installation, see documentation of pip.'
        test_class.fail(msg)
#    test_class.assertTrue(protobuf_installed)


def optimize_in_chunks( sd, sim_time_ns, ns_per_chunk ):
    """
        Optimizes sd->bd() in nchunks iterations, writing statistics at
        the end of each iteration

        returns the rest lengths of the all chain bonds at the end of each chunk,
        for bonds of RelaxingSpring type, based on their RMF values
        """
    dT_fs= sd.get_bd().get_maximum_time_step()
    fs_per_chunk= ns_per_chunk*1E+6
    nframes_per_chunk= int(round(fs_per_chunk/dT_fs))
    nframes_per_chunk= max(nframes_per_chunk, 1)
    sim_time_fs= sim_time_ns*1E+6
    nframes= int(round(sim_time_fs/dT_fs))
    nframes= max(nframes, 1)
    timer = IMP.npctransport.create_boost_timer()
    nframes_left = nframes
    R=[]
    while(nframes_left > 0):
        nframes_per_chunk = min(nframes_per_chunk, nframes_left)
        sd.get_bd().optimize( nframes_per_chunk )
        ps=sd.get_beads()
        cur_R=[]
        for p in ps:
            if IMP.npctransport.RelaxingSpring.get_is_setup(p):
                rs= IMP.npctransport.RelaxingSpring(p)
                cur_R.append(rs.get_rest_length())
        R.append(cur_R)
        if(sd.get_statistics().get_is_activated()):
            sd.get_statistics().update( timer, nframes - nframes_per_chunk ) # TODO: timer?
        nframes_left = nframes_left - nframes_per_chunk
    return np.array(R)

def check_import_pandas_with_series_autocorr():
    '''
    import pandas, return true if successful and pandas.Series
    has an autocorr() method, return false otherwise
    '''
    try:
        import pandas as pd
    except ImportError:
        return False
    autocorr_func = getattr(pd.Series, "autocorr", None)
    return callable(autocorr_func)
