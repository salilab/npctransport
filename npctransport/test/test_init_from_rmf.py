import IMP
import IMP.test
import RMF
import IMP.rmf
import IMP.container
import math
import IMP.base
from IMP.npctransport import *

class Tests(IMP.test.TestCase):

    def _make_and_run_simulation(self, output_pb_file, output_rmf_file, seed):
        IMP.base.random_number_generator.seed( seed )
        config= IMP.npctransport.get_data_path( "quick2.pb" );
        IMP.set_log_level( IMP.SILENT );
        print "assigning parameter ranges"
        num=assign_ranges( config, output_pb_file,
                          0, True, seed );
        print "num ranges %d" % num
        sd = SimulationData( output_pb_file, False, output_rmf_file );
        sd.get_m().set_log_level( IMP.SILENT );
        print "Files are " + output_pb_file + \
              " and " + output_rmf_file
        timer = IMP.npctransport.create_boost_timer()
        sd.get_bd().optimize( 1000 )
        # make sure final state is written
        sd.get_rmf_sos_writer().update_always()
        # make sure statistics with final rmf conformation is written
        sd.update_statistics( timer )
        return sd

    def _get_diffuser_coords( self, sd, i ):
        """ returns the coordinates of diffuser number i in SimulationData
            object sd """
        d = sd.get_diffusers().get_particles()[i]
        d_coords = IMP.core.XYZ( d ).get_coordinates()
        d_rframe = IMP.core.RigidBody( d ).get_reference_frame()
        print "Diffuser refframe" + str(i) + ":", d,  d_rframe
        return d_coords

    def _get_diffusers_coords( self, sd ):
        """ get the coordinates list of all diffusers in SimulationData
            object sd """
        coords = []
        for i in range( sd.get_diffusers().get_number_of_particles() ):
            coords.append( self._get_diffuser_coords( sd, i ) )
        return coords

    def test_init_from_rmf(self):
        """ Testing whether initialize_positions_from_rmf indeed
            restores the diffusers coordinates correctly """
        # First simulation with one seed:
        output_pb1= self.get_tmp_file_name( "output1.pb" );
        output_rmf1= IMP.base.create_temporary_file_name( "output", ".rmf" );
        print "Starting first simulation with RMF file " + output_rmf1
        sd = self._make_and_run_simulation( output_pb1, output_rmf1, seed = 1)
        e1 = sd.get_bd().get_scoring_function().evaluate(False)
        print "*** After first simulation"
        print "Energy %.2f" % e1
        coordsI = self._get_diffusers_coords( sd )
        sd = None

        # Second simulation with another seed:
        output_pb2= self.get_tmp_file_name( "output2.pb" );
        output_rmf2= IMP.base.create_temporary_file_name( "output2", ".rmf" );
        print "*** Starting second simulation with RMF file " + output_rmf2
        sd = self._make_and_run_simulation( output_pb2, output_rmf2, seed = 2)
        e2 = sd.get_bd().get_scoring_function().evaluate(False)
        print "*** After second simulation"
        print "Energy %.2f" % e2
        coordsII = self._get_diffusers_coords( sd )

        # Restoration from first simulation through rmf file:
        print "*** Initializing positions from RMF file " + output_rmf1
        fl= RMF.open_rmf_file_read_only(output_rmf1)
        sd.initialize_positions_from_rmf( fl )
        e3 = sd.get_bd().get_scoring_function().evaluate(False)
        print "*** After initializing positions from RMF file " + output_rmf1
        print "Energy %.2f" % e3
        self.assertEqual(e1, e3)
        self.assertNotEqual(e1, e2)
        coordsIII = self._get_diffusers_coords( sd )
        # make sure that all coordinates were restored and as sanity
        # check control, also that the second simulation is different than
        # the first one:
        for i in range( len( coordsI ) ):
            for j in range(3):
                self.assertEqual(coordsI[i][j], coordsIII[i][j])
                self.assertNotEqual(coordsI[i][j], coordsII[i][j])

        # Simulate more to scramble stuff
        sd.get_bd().optimize( 1000 )
        e4 = sd.get_bd().get_scoring_function().evaluate(False)
        print "Energy %.2f" % e4
        coordsIV = self._get_diffusers_coords( sd )

        # Restoration from first simulation through output protobuf file:
        print "*** Initializing positions from ProtoBuf file " + output_pb1
        f=open(output_pb1, "rb")
        config= Output()
        config.ParseFromString(f.read())
        fl= RMF.open_rmf_buffer_read_only( config.rmf_conformation )
        sd.initialize_positions_from_rmf( fl )
        e5 = sd.get_bd().get_scoring_function().evaluate(False)
        print "*** After initializing positions from RMF file " + output_rmf1
        print "Energy %.2f" % e5
        self.assertEqual(e1, e5)
        self.assertNotEqual(e1, e4)
        coordsV = self._get_diffusers_coords( sd )
        for i in range( len( coordsI ) ):
            for j in range(3):
                self.assertEqual(coordsI[i][j], coordsV[i][j])
                self.assertNotEqual(coordsI[i][j], coordsIV[i][j])


if __name__ == '__main__':
    IMP.test.main()
