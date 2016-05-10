from __future__ import print_function
import IMP
import IMP.test
import RMF
import IMP.rmf
import IMP.container
import math
import IMP
import test_util
from IMP.npctransport import *
import time

# TODO: currently this is slow only due to initial optimization algorithm
#       that must be imporved


class Tests(IMP.test.TestCase):

    def run_from_config(self, config, output):
        """
        Run using work-unit 0 using specified config file,
        dumping output to specified output file

        return - the resulting simulation data object
        """
        IMP.set_log_level(IMP.SILENT)
        print("assigning parameter ranges from config")
        num = assign_ranges(config, output,
                            0, True, 10)
#        IMP.set_log_level(IMP.TERSE)
        sd = IMP.npctransport.SimulationData(output, False)
        sd.set_rmf_file( self.get_tmp_file_name("out0.rmf"), False )
        print("BEFORE INIT", time.ctime())
        if IMP.get_check_level() >= IMP.USAGE_AND_INTERNAL:
            short_init_factor = 0.00001
            opt_cycles = 2
        else:
            short_init_factor = 0.01
            opt_cycles = 10000
        sd.get_bd().set_log_level(IMP.SILENT)
        IMP.npctransport.initialize_positions(sd, [], False, short_init_factor)
        print("AFTER INIT", time.ctime())
        sd.activate_statistics()
        obd = sd.get_bd()
        obd.optimize(opt_cycles)
        print("AFTER OPTIMIZE", time.ctime())
        timer = IMP.npctransport.timer()
        # lame test
        # rt= sd.get_root()
        # rtt= IMP.npctransport.Transporting.setup_particle(rt, True)
        # rtf= rt.get_child(0)
        # rttf= IMP.npctransport.Transporting.setup_particle(rtf, False)
        print("updating stats")
        sd.get_statistics().update(timer, 0)
        return sd

    def assert_transporting_equal(self, sd1, sd2):
        """ assert that sd1 and sd2 have identical Transporting statistics """
        for d1, d2 in zip(sd1.get_beads(),
                          sd2.get_beads()):
            if(not IMP.npctransport.Transporting.get_is_setup(d1)):
                continue
            if(not IMP.npctransport.Transporting.get_is_setup(d2)):
                continue
            t1 = IMP.npctransport.Transporting(d1)
            t2 = IMP.npctransport.Transporting(d2)
            print ("Bead particles: ")
            print (d1, d2)
            print ("Comparing transport statistics: ", t1, t2)
            self.assert_(t1.get_is_last_entry_from_top()
                         == t2.get_is_last_entry_from_top())
            self.assertAlmostEqual(t1.get_last_tracked_z(),
                                   t2.get_last_tracked_z(), delta=.00001)
            self.assert_(t1.get_n_entries_bottom()
                         == t2.get_n_entries_bottom())
            self.assert_(t1.get_n_entries_top()
                         == t2.get_n_entries_top())

    def assert_almost_equal_sds(self, sd1, sd2):
        """
        assert that sd1 and sd2 has nearly identical positions for beads
        and sites + identical timers and Transporting porperties
        """
        # check beads refframes
        for p, pp in zip(sd1.get_beads(),
                         sd2.get_beads()):
            self.assert_((IMP.core.XYZ(p).get_coordinates()
                          - IMP.core.XYZ(pp).get_coordinates()).get_magnitude() < .0001)
            q0 = IMP.core.RigidBody(
                p).get_reference_frame(
            ).get_transformation_to(
            ).get_rotation(
            ).get_quaternion(
            )
            q1 = IMP.core.RigidBody(
                pp).get_reference_frame(
            ).get_transformation_to(
            ).get_rotation(
            ).get_quaternion(
            )
            print(q0, q1)
            for qa, qb in zip(q0, q1):
                self.assertAlmostEqual(qa, qb, delta=.01)
        # check sites
        sites0 = sd1.get_sites(IMP.core.ParticleType("kap"))
        sites1 = sd2.get_sites(IMP.core.ParticleType("kap"))
        for s0, s1 in zip(sites0, sites1):
            self.assert_(IMP.algebra.get_distance(s0, s1) < .0001)
        # check timers
        bd1 = sd1.get_bd()
        bd2 = sd2.get_bd()
        self.assert_(bd2.get_current_time() > 0)
        self.assert_(bd1.get_current_time() == bd2.get_current_time())
        # check Transporting
        self.assert_transporting_equal(sd1, sd2)

    def test_init_from_output(self):
        """ Testing whether positions are loaded properly from output file """
        print("TEST_INIT_FROM_OUTPUT")
        # random generator initialization
        IMP.set_log_level(IMP.SILENT)
        config = self.get_tmp_file_name("simple_cfg.pb")
        test_util.make_simple_cfg(
            config,
            is_slab_on=True,
            n_particles_factor=1.5)
        rt_output = self.get_tmp_file_name("round_trip_output.pb")
        print("RT output: ", rt_output)
        sd = self.run_from_config(config, rt_output)

        print("reloading from output file ", rt_output)
        sdp = IMP.npctransport.SimulationData(rt_output, False)
        sdp.activate_statistics()
        sd.set_rmf_file( self.get_tmp_file_name("out1.rmf"), False )
        print("After reload", time.ctime())
        self.assert_almost_equal_sds(sd, sdp)
#        print "updating stats at end"
#        sd.update_statistics(timer, 0);

    def test_init_from_old_output1(self):
        """ Testing whether an old output file is loaded properly """
        print("TEST_INIT_FROM_OLD_OUTPUT1")
        expected_sites = [ (3.67394e-15, 0, -30),
                           (17.2447, -0.377296, -24.5455),
                           (-1.55764, -23.0892, -19.0909),
                           (-19.7675, 17.9804, -13.6364),
                           (28.589, -3.96571, -8.18182),
                           (-28.2097, -9.8375, -2.72727),
                           (20.6227, 21.6163, 2.72727),
                           (-4.64461, -28.4866, 8.18182),
                           (-17.9874, 19.7612, 13.6364),
                           (17.8043, 14.7833, 19.0909),
                           (13.5084, 10.7259, 24.5455),
                           (0, 0, 30) ]
        expected_time =  12500039289
        expected_particles=[ [-253.636, -108.652, 40.4134, 1, 0, 0, 0], #trans + quaternion
                             [-236.08, -127.91, 27.47, 0.23, 0.33, -0.84, -0.36] ]

        # random generator initialization
        # RMF.set_log_level("trace")
        IMP.set_log_level(IMP.SILENT)
        rt_prev_output = self.get_input_file_name("out149.pb")
        out_rmf = self.get_tmp_file_name("movie.rmf")
        rt_new_output = self.get_tmp_file_name("out.pb")
        print("reloading from output file ", rt_prev_output)
        print("New Output: ", rt_new_output)
        sd = IMP.npctransport.SimulationData(rt_prev_output,
                                             False,
                                             out_rmf,
                                             rt_new_output)
        sd.activate_statistics()
        sd.set_rmf_file(out_rmf, False)
        for i,p in enumerate(sd.get_beads()):
            if( i >= len(expected_particles) ):
                break
            pc= [x for x in IMP.core.XYZ(p).get_coordinates()]
            pq = [x  for x in IMP.core.RigidBody(p
               ).get_reference_frame(
               ).get_transformation_to(
               ).get_rotation(
               ).get_quaternion() ]
            for x,y in zip(pc+pq, expected_particles[i]):
                self.assertAlmostEqual(x,y, delta=0.1)
#            if(not IMP.npctransport.Transporting.get_is_setup(p)):
#                continue
#            pt = IMP.npctransport.Transporting(p)
#            print "T: ", pt.get_is_last_entry_from_top(),
#            print pt.get_last_tracked_z(),
#            print pt.get_n_entries_bottom(),
#            print pt.get_n_entries_top(),
#            print
        sites = sd.get_sites(IMP.core.ParticleType("kap"))
        for i,s in enumerate(sites):
            self.assertAlmostEqual(
                ( s.get_center() - expected_sites[i]).get_magnitude(),
                0, delta=.01 )
        # check timers
        t = sd.get_bd().get_current_time()
        self.assertAlmostEqual(t, expected_time, delta=1)

    def test_init_from_old_output_more_recent(self):
        """
        Testing whether a more recent output file (Dec 2013) is loaded
        at all
        """
        print("TEST_INIT_FROM_OLD_OUTPUT_MORE_RECENT")
        IMP.set_log_level(IMP.SILENT)
        rt_prev_output = self.get_input_file_name("out_more_recent.pb")
        out_rmf = self.get_tmp_file_name("movie.rmf")
        rt_new_output = self.get_tmp_file_name("out.pb")
        print("reloading from output file ", rt_prev_output)
        print("New Output: ", rt_new_output)
        sd = IMP.npctransport.SimulationData(rt_prev_output,
                                             False,
                                             out_rmf,
                                             rt_new_output)
        sd.activate_statistics()
        sd.set_rmf_file(out_rmf, False)
#        sd.get_bd().optimize(1)




if __name__ == '__main__':
    IMP.test.main()
