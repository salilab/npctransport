import IMP
import IMP.test
import RMF
import IMP.rmf
import IMP.container
import math
import IMP.base
from IMP.npctransport import *

class Tests(IMP.test.TestCase):

    def _make_and_run_simulation(self, output_file, seed):
        RMF.set_show_hdf5_errors( True );
        config= IMP.npctransport.get_data_path( "quick.pb" );
        assignment= self.get_tmp_file_name( "output.pb" );
        IMP.set_log_level( IMP.SILENT );
        print "assigning parameter ranges"
        num=assign_ranges( config, assignment,
                          0, True, seed );
        print "num ranges %d" % num
        sd = SimulationData( assignment, False, output_file );
        sd.get_m().set_log_level( IMP.SILENT );
        print "Files are " + assignment + \
              " and " + output_file
        sd.get_bd().optimize( 10 )
        # make sure final state is written
        sd.get_rmf_sos_writer().update_always()
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
        # random generator initialization
        seed = 1 # use time instead?
        IMP.base.random_number_generator.seed( seed )

        # First simulation:
        output= IMP.base.create_temporary_file_name( "output", ".rmf" );
        print "Starting first simulation with RMF file " + output
        sd = self._make_and_run_simulation( output, seed )
        print "*** After first simulation"
        coordsI = self._get_diffusers_coords( sd )

        # Second simulation:
        output2= IMP.base.create_temporary_file_name( "output2", ".rmf" );
        print "*** Starting second simulation with RMF file " + output2
        sd2 = self._make_and_run_simulation( output2, seed )
        print "*** After second simulation"
        coordsII = self._get_diffusers_coords( sd2 )

        # Restoration from first simulation:
        print "*** Initializing positions from RMF file " + output
        fl= RMF.open_rmf_file_read_only(output)
        sd2.initialize_positions_from_rmf( fl )
        print "*** After initializing positions from RMF file " + output
        coordsIII = self._get_diffusers_coords( sd2 )

        # make sure that all coordinates were restored and as sanity
        # check control, also that the second simulation is different than
        # the first one:
        for i in range( len( coordsI ) ):
            for j in range(3):
                print i, j
                self.assertEqual(coordsI[i][j], coordsIII[i][j])
                self.assertNotEqual(coordsI[i][j], coordsII[i][j])


if __name__ == '__main__':
    IMP.test.main()
