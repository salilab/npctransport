from __future__ import print_function
import IMP
import IMP.test
import IMP.algebra
import IMP.npctransport
import IMP.container
import math
import test_util


def create_reference_frame(trans_vec, rot_axis, rot_angle_radians):
    translation= [5,-5,7]
    rotation= IMP.algebra.get_rotation_about_axis(rot_axis, rot_angle_radians)
    transformation= IMP.algebra.Transformation3D(rotation, trans_vec)
    ref_frame= IMP.algebra.ReferenceFrame3D(transformation)
    return ref_frame

def are_rigid_body_reference_frames_equal(rb1, rb2):
    EPSILON= 1E-10
    T1= rb1.get_reference_frame().get_transformation_from()
    T2= rb2.get_reference_frame().get_transformation_to()
    T= T1*T2
    # If equal, T should be identity, so check all base vectors
    # are transformed to themselves
    for v in ([1,0,0],[0,1,0],[0,0,1]):
        v3= IMP.algebra.Vector3D(v)
        if (T*v3-v3).get_magnitude() > EPSILON:
            return False
    return True

def get_fg_beads_from_simulation_data(sd):
    ret_value=[]
    for p in sd.get_beads():
        pi= p.get_index()
        if sd.get_is_fg_bead(pi):
            ret_value.append(p)
    return ret_value


def run_from_config(config_fname, output_fname, rmf_fname=None):
        """
        Run using work-unit 0 using specified config file,
        dumping output to specified output file

        return - the resulting simulation data object
        """
        IMP.set_log_level(IMP.SILENT)
        print("assigning parameter ranges from config into output", output_fname)
        num = IMP.npctransport.assign_ranges(config_fname, output_fname,
                            0, True, 10)
        sd = IMP.npctransport.SimulationData(output_fname, False)
        sd.set_log_level(IMP.SILENT)
        if rmf_fname is not None:
            sd.set_rmf_file( rmf_fname, False )
        if IMP.get_check_level() >= IMP.USAGE_AND_INTERNAL:
            short_init_factor = 0.00001
            opt_cycles = 1
        else:
            short_init_factor = 0.01
            opt_cycles = 10000
        sd.get_bd().set_log_level(IMP.SILENT)
        IMP.npctransport.initialize_positions(sd, [], False, short_init_factor)
        sd.activate_statistics()
        obd = sd.get_bd()
        obd.optimize(opt_cycles)
        timer = IMP.npctransport.timer()
        # lame test
        # rt= sd.get_root()
        # rtt= IMP.npctransport.Transporting.setup_particle(rt, True)
        # rtf= rt.get_child(0)
        # rttf= IMP.npctransport.Transporting.setup_particle(rtf, False)
        print("updating stats")
        sd.get_statistics().update(timer, 0)
        return sd


class FGCopyTests(IMP.test.TestCase):
    def test_copy_particle_reference_frame(self):
        ''' test utility function for copying one particle's xyz coordinates and reference frame
            to another
        '''
        # XYZ:
        m= IMP.Model()
        src_p= IMP.Particle(m)
        trg_p= IMP.Particle(m)
        src_xyz= IMP.core.XYZ.setup_particle(src_p, [3,1,1])
        trg_xyz= IMP.core.XYZ.setup_particle(trg_p, [0,0,0])
        self.assertNotAlmostEqual(src_xyz.get_vector_to(trg_xyz).get_magnitude(),
                                  0.0)
        IMP.npctransport.copy_particle_reference_frame_if_applicable(src_p, trg_p)
        self.assertAlmostEqual(src_xyz.get_vector_to(trg_xyz).get_magnitude(),
                               0.0)
        # RigidBody:
        rf1= create_reference_frame([1,1,2],[0,1,1], math.pi/2.0)
        rf2= create_reference_frame([-1,-1,3],[1,3,-1], math.pi/4.0)
        src_rb= IMP.core.RigidBody.setup_particle(src_p, rf1)
        if IMP.get_check_level() >= IMP.USAGE:
            with self.assertRaises(IMP.UsageException):
                IMP.npctransport.copy_particle_reference_frame_if_applicable(src_p,
                                                                             trg_p)
        trg_rb= IMP.core.RigidBody.setup_particle( trg_p, rf2)
        self.assertFalse(are_rigid_body_reference_frames_equal(src_rb, trg_rb))
        IMP.npctransport.copy_particle_reference_frame_if_applicable(src_p, trg_p)
        self.assertTrue(are_rigid_body_reference_frames_equal(src_rb, trg_rb))


    def test_copy_fgs(self):
        test_util.test_protobuf_installed(self)
        config_fname = "./config_fg_copy.pb" # self.get_tmp_file_name("simple_cfg.pb")
        output_fname1 = "./output_fg_copy.pb" # self.get_tmp_file_name("test_fg_copy.pb")
        output_fname2 = "./output_fg_copy.pb" # self.get_tmp_file_name("test_fg_copy.pb")
        rmf_fname1 = "./rmf_fg_copy.rmf" # self.get_tmp_file_name("test_fg_copy.pb")
        rmf_fname2 = "./rmf_fg_copy2.rmf" # self.get_tmp_file_name("test_fg_copy.pb")
        test_util.make_simple_cfg(
            config_fname,
            is_slab_on=True,
            n_particles_factor=3.0,
            is_obstacles=True)
        sd1= run_from_config(config_fname, output_fname1, rmf_fname1)
        sd2= run_from_config(config_fname, output_fname2, rmf_fname2)
        fg_beads1= get_fg_beads_from_simulation_data(sd1)
        fg_beads2= get_fg_beads_from_simulation_data(sd2)
        # Verify that not all are equal initially:
        are_all_equal= True
        for p1, p2 in zip(fg_beads1, fg_beads2):
            if not are_rigid_body_reference_frames_equal \
               ( IMP.core.RigidBody(p1), IMP.core.RigidBody(p2) ):
                are_all_equal=False
        self.assertFalse(are_all_equal)
        # Copy then verify that all are equal now:
        IMP.npctransport.copy_FGs_coordinates(sd1, sd2)
        for p1, p2 in zip(fg_beads1, fg_beads2):
            self.assertTrue(are_rigid_body_reference_frames_equal \
                            ( IMP.core.RigidBody(p1), IMP.core.RigidBody(p2) ) )
        if IMP.get_check_level() >= IMP.USAGE_AND_INTERNAL:
            sd2.get_bd().optimize(1)
        else:
            sd2.get_bd().optimize(20000)



if __name__ == '__main__':
    IMP.test.main()
