import IMP
import IMP.test
import RMF
import IMP.rmf
import IMP.container
import math
import IMP.base
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
        IMP.set_log_level( IMP.SILENT );
        print "assigning parameter ranges from config"
        num=assign_ranges( config, output,
                          0, True, 10 );
#        IMP.set_log_level(IMP.TERSE)
        sd= IMP.npctransport.SimulationData(output, False,
                                            self.get_tmp_file_name("out0.rmf"));
        print "BEFORE INIT", time.ctime()
        if IMP.base.get_check_level() >= IMP.base.USAGE_AND_INTERNAL:
            short_init_factor=0.00001
            opt_cycles = 2
        else:
            short_init_factor=0.01
            opt_cycles = 10000
        IMP.npctransport.initialize_positions(sd, [], False, short_init_factor)
        print "AFTER INIT", time.ctime()
        obd= sd.get_bd()
        obd.optimize(opt_cycles)
        print "AFTER OPTIMIZE", time.ctime()
        timer= IMP.npctransport.timer();
        # # lame test
        # rt= sd.get_root()
        # rtt= IMP.npctransport.Transporting.setup_particle(rt, True)
        # rtf= rt.get_child(0)
        # rttf= IMP.npctransport.Transporting.setup_particle(rtf, False)
        print "updating stats"
        sd.get_statistics().update(timer, 0);
        return sd

    def assert_transporting_equal(self, sd1, sd2):
        """ assert that sd1 and sd2 have identical Transporting statistics """
        for d1, d2 in zip(sd1.get_diffusers().get_particles(),
                         sd2.get_diffusers().get_particles()):
           if( not IMP.npctransport.Transporting.particle_is_instance(d1) ):
               continue
           if( not IMP.npctransport.Transporting.particle_is_instance(d2) ):
               continue
           t1 = IMP.npctransport.Transporting( d1 )
           t2 = IMP.npctransport.Transporting( d2 )
           print "Diffuser particles: "
           print d1, d2
           print "Comparing transport statistics: ", t1, t2
           self.assert_(t1.get_is_last_entry_from_top()
                        == t2.get_is_last_entry_from_top() )
           self.assertAlmostEqual(t1.get_last_tracked_z(),
                        t2.get_last_tracked_z(), delta=.00001 )
           self.assert_(t1.get_n_entries_bottom()
                        == t2.get_n_entries_bottom() );
           self.assert_(t1.get_n_entries_top()
                        == t2.get_n_entries_top() );

    def assert_almost_equal_sds(self, sd1, sd2):
        """
        assert that sd1 and sd2 has nearly identical positions for diffusers
        and sites + identical timers and Transporting porperties
        """
        # check diffusers refframes
        for p, pp in zip(sd1.get_diffusers().get_particles(),
                         sd2.get_diffusers().get_particles()):
            self.assert_((IMP.core.XYZ(p).get_coordinates()
                          - IMP.core.XYZ(pp).get_coordinates()).get_magnitude() < .0001)
            q0= IMP.core.RigidBody(p).get_reference_frame().get_transformation_to().get_rotation().get_quaternion()
            q1= IMP.core.RigidBody(pp).get_reference_frame().get_transformation_to().get_rotation().get_quaternion()
            print q0, q1
            for qa, qb in zip(q0, q1):
                self.assertAlmostEqual(qa, qb, delta=.01)
        # check sites
        sites0= sd1.get_sites(IMP.core.ParticleType("kap"))
        sites1= sd2.get_sites(IMP.core.ParticleType("kap"))
        for s0,s1 in zip(sites0, sites1):
            self.assert_(IMP.algebra.get_distance(s0,s1) < .0001)
        # check timers
        bd1 = sd1.get_bd()
        bd2 = sd2.get_bd()
        self.assert_(bd2.get_current_time() >0)
        self.assert_(bd1.get_current_time()==bd2.get_current_time())
        # check Transporting
        self.assert_transporting_equal(sd1, sd2)

    def test_init_from_output(self):
        """ Testing whether positions are loaded properly from output file """
        # random generator initialization
        IMP.set_log_level(IMP.SILENT)
        config= self.get_tmp_file_name( "simple_cfg.pb") ;
        test_util.make_simple_cfg( config, is_slab_on = True, n_particles_factor = 1.5)
        rt_output= self.get_tmp_file_name("round_trip_output.pb")
        print "RT output: ", rt_output
        sd = self.run_from_config( config, rt_output )

        print "reloading from output file ", rt_output
        sdp= IMP.npctransport.SimulationData(rt_output, False,
                                             self.get_tmp_file_name("out1.rmf"));
        print "After reload", time.ctime()
        self.assert_almost_equal_sds(sd, sdp)
#        print "updating stats at end"
#        sd.update_statistics(timer, 0);

if __name__ == '__main__':
    IMP.test.main()
