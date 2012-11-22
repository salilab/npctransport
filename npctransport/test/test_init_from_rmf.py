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
        RMF.set_show_hdf5_errors( True );
        config= IMP.npctransport.get_data_path( "quick.pb" );
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
        sd.get_bd().optimize( 10 )
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
        print "Diffuser " + str(i) + ":", d, d_coords
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
        print "*** After first simulation"
        coordsI = self._get_diffusers_coords( sd )

        # Second simulation with another seed:
        output_pb2= self.get_tmp_file_name( "output2.pb" );
        output_rmf2= IMP.base.create_temporary_file_name( "output", ".rmf" );
        print "*** Starting second simulation with RMF file " + output_rmf2
        sd2 = self._make_and_run_simulation( output_pb2, output_rmf2, seed = 2)
        print "*** After second simulation"
        coordsII = self._get_diffusers_coords( sd2 )

        # Restoration from first simulation through rmf file:
        print "*** Initializing positions from RMF file " + output_rmf1
        fl= RMF.open_rmf_file_read_only(output_rmf1)
        sd2.initialize_positions_from_rmf( fl )
        print "*** After initializing positions from RMF file " + output_rmf1
        coordsIII = self._get_diffusers_coords( sd2 )
        # make sure that all coordinates were restored and as sanity
        # check control, also that the second simulation is different than
        # the first one:
        for i in range( len( coordsI ) ):
            for j in range(3):
                self.assertEqual(coordsI[i][j], coordsIII[i][j])
                self.assertNotEqual(coordsI[i][j], coordsII[i][j])

        # Simulate more to scramble stuff
        sd2.get_bd().optimize( 20 )
        coordsIV = self._get_diffusers_coords( sd2 )

        # Restoration from first simulation through output protobuf file:
        print "*** Initializing positions from ProtoBuf file " + output_pb1
        f=open(output_pb1, "rb")
        config= Output()
        config.ParseFromString(f.read())
        fl= RMF.open_rmf_buffer_read_only( config.rmf_conformation )
        sd.initialize_positions_from_rmf( fl )
        print "*** After initializing positions from RMF file " + output_rmf1
        coordsV = self._get_diffusers_coords( sd )
        for i in range( len( coordsI ) ):
            for j in range(3):
                self.assertEqual(coordsI[i][j], coordsV[i][j])
                self.assertNotEqual(coordsI[i][j], coordsIV[i][j])


if __name__ == '__main__':
    IMP.test.main()
