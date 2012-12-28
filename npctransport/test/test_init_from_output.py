import IMP
import IMP.test
import RMF
import IMP.rmf
import IMP.container
import math
import IMP.base
from IMP.npctransport import *

class Tests(IMP.test.TestCase):
    def test_init_from_rmf(self):
        """ Testing whether positions are loaded properly from output file """
        # random generator initialization
        IMP.set_log_level(IMP.SILENT)
        config= IMP.npctransport.get_data_path( "quick.pb" );
        output= self.get_tmp_file_name("round_trip_output.pb")
        IMP.set_log_level( IMP.SILENT );
        print "assigning parameter ranges"
        num=assign_ranges( config, output,
                          0, True, 10 );
        sd= IMP.npctransport.SimulationData(output, False,
                                            self.get_tmp_file_name("out0.rmf"));
        IMP.npctransport.initialize_positions(sd, [], False)
        obd= sd.get_bd()
        obd.optimize(1)
        timer= IMP.npctransport.timer();
        # lame test
        rt= sd.get_root()
        rtt= IMP.npctransport.Transporting.setup_particle(rt, True)
        rtf= rt.get_child(0)
        rttf= IMP.npctransport.Transporting.setup_particle(rtf, False)
        print "updating stats"
        sd.update_statistics(timer, 0);
        sites0= sd.get_sites(IMP.core.ParticleType("kap"))

        print "reloading"
        sdp= IMP.npctransport.SimulationData(output, False,
                                             self.get_tmp_file_name("out1.rmf"));
        sites1= sd.get_sites(IMP.core.ParticleType("kap"))
        for p, pp in zip(sd.get_diffusers().get_particles(),
                         sdp.get_diffusers().get_particles()):
            self.assert_((IMP.core.XYZ(p).get_coordinates()
                         - IMP.core.XYZ(pp).get_coordinates()).get_magnitude() < .0001)
            q0= IMP.core.RigidBody(p).get_reference_frame().get_transformation_to().get_rotation().get_quaternion()
            q1= IMP.core.RigidBody(pp).get_reference_frame().get_transformation_to().get_rotation().get_quaternion()
            print q0, q1
            for qa, qb in zip(q0, q1):
                self.assertAlmostEqual(qa, qb, delta=.01)
        for s0,s1 in zip(sites0, sites1):
            self.assert_(IMP.algebra.get_distance(s0,s1) < .0001)
        bd = sdp.get_bd()
        self.assert_(bd.get_current_time() >0)
        self.assert_(bd.get_current_time()==obd.get_current_time())
        rt= sdp.get_root()
        self.assert_(rtt.get_is_last_entry_from_top())
        rtf= rt.get_child(0)
        self.assert_(not rttf.get_is_last_entry_from_top());
        print "updating stats at end"
        sd.update_statistics(timer, 0);

if __name__ == '__main__':
    IMP.test.main()
