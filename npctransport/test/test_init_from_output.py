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
        sites0= sd.get_sites(IMP.core.ParticleType("kap"))
        IMP.npctransport.initialize_positions(sd, [], False)
        timer= IMP.npctransport.timer();
        print "updating stats"
        sd.update_statistics(timer, 0);

        print "reloading"
        sdp= IMP.npctransport.SimulationData(output, False,
                                             self.get_tmp_file_name("out1.rmf"));
        sites1= sd.get_sites(IMP.core.ParticleType("kap"))
        for p, pp in zip(sd.get_diffusers().get_particles(),
                         sdp.get_diffusers().get_particles()):
            self.assert_((IMP.core.XYZ(p).get_coordinates()
                         - IMP.core.XYZ(pp).get_coordinates()).get_magnitude() < .0001)
        for s0,s1 in zip(sites0, sites1):
            self.assert_(IMP.algebra.get_distance(s0,s1) < .0001)


if __name__ == '__main__':
    IMP.test.main()
