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
        config= IMP.npctransport.get_data_path( "quick.pb" );
        assignment= self.get_tmp_file_name( "assignment.pb" );
        IMP.set_log_level( IMP.SILENT );
        print "assigning parameter ranges"
        num=assign_ranges( config, assignment,
                          0, True, seed );
        output= self.get_tmp_file_name("round_trip_output.pb")
        sd= IMP.npctranspot.SimulationData(assignment, False, output);
        IMP.npctransport.initialize_positions(sd, [], False)
        timer= IMP.npctransport.timer();
        sd.update_statistics(timer, 0);

        sdp= IMP.npctranspot.SimulationData(output, False, output);

        for p, pp in zip(sd.get_diffusers().get_particles(),
                         sdp.get_diffusers().get_particles()):
            self.assert_((IMP.core.XYZ(p).get_coordinates()
                         - IMP.core.XYZ(pp).get_coordinates()).get_magnitude() < .0001)


if __name__ == '__main__':
    IMP.test.main()
