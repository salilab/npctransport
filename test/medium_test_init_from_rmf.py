from __future__ import print_function
import IMP
import IMP.test
import RMF
import IMP.rmf
import IMP.container
import math
import IMP
from IMP.npctransport import *
import test_util

class Tests(IMP.test.TestCase):

    def _make_and_run_simulation(self, output_pb_file, output_rmf_file, seed):
        IMP.random_number_generator.seed( seed )
        config= IMP.npctransport.get_data_path( "quick2.pb" )
        IMP.set_log_level( IMP.SILENT );
        print("assigning parameter ranges")
        num=assign_ranges( config, output_pb_file,
                          0, True, seed );
        print("num ranges %d" % num)
        sd = SimulationData( output_pb_file, False )
        sd.activate_statistics()
        sd.set_rmf_file(output_rmf_file, False)
        sd.get_model().set_log_level( IMP.SILENT );
        print("Files are " + output_pb_file + \
              " and " + output_rmf_file)
        timer = IMP.npctransport.create_boost_timer()
        sd.get_bd().optimize( 1000 )
        # make sure final state is written to rmf
        sd.get_rmf_sos_writer().update_always()
        # write statistics and final rmf conformation
        sd.get_statistics().update( timer )
        return sd

    def _get_bead_coords( self, p ):
        """ returns the coordinates of the particle """
        rb = IMP.core.RigidBody( p )
        coords = rb.get_coordinates()
        refframe = rb.get_reference_frame()
        print("Bead " + str(p) + "," + str(p.get_index()) + ":",   refframe)
        return coords

    def _get_beads_coords( self, sd ):
        """ get the coordinates list of all beads in SimulationData
            object sd """
        coords = []
        for bead in sd.get_beads():
            coords.append( self._get_bead_coords( bead ) )
        return coords

    def test_init_from_rmf(self):
        """ Testing whether initialize_positions_from_rmf indeed
            restores the beads coordinates correctly """
        IMP.set_log_level( IMP.SILENT );
        test_util.test_protobuf_installed(self)
        # First simulation with one seed:
        output_pb1= self.get_tmp_file_name( "output1.pb" );
        output_rmf1= IMP.create_temporary_file_name( "output", ".rmf" );
        print("Starting first simulation with RMF file " + output_rmf1)
        sd = self._make_and_run_simulation( output_pb1, output_rmf1, seed = 1)
        e1 = sd.get_bd().get_scoring_function().evaluate(False)
        print ("*** After first simulation")
        print ("Energy %.2f" % e1)
        coordsI = self._get_beads_coords( sd )
        sd = None

        # Second simulation with another seed:
        output_pb2= self.get_tmp_file_name( "output2.pb" );
        output_rmf2= IMP.create_temporary_file_name( "output2", ".rmf" );
        print("*** Starting second simulation with RMF file " + output_rmf2)
        print("YYY")
        sd = self._make_and_run_simulation( output_pb2, output_rmf2, seed = 2)
        print("XXX")
        e2 = sd.get_bd().get_scoring_function().evaluate(False)
        print ("*** After second simulation")
        print ("Energy %.2f" % e2)
        coordsII = self._get_beads_coords( sd )

        # Restoration from first simulation through rmf file:
        print("*** Initializing positions from RMF file " + output_rmf1)
        fl= RMF.open_rmf_file_read_only(output_rmf1)
        sd.initialize_positions_from_rmf( fl )
        e3 = sd.get_bd().get_scoring_function().evaluate(False)
        print("*** After initializing positions from RMF file " + output_rmf1)
        print("Energy %.2f" % e3)
        self.assertAlmostEqual(e1, e3, delta=0.0001)
        self.assertNotAlmostEqual(e1, e2, delta=0.0001)
        coordsIII = self._get_beads_coords( sd )
        # make sure that all coordinates were restored and as sanity
        # check control, also that the second simulation is different than
        # the first one:
        for i in range( len( coordsI ) ):
            for j in range(3):
                self.assertAlmostEqual(coordsI[i][j], coordsIII[i][j], delta=0.00001)
                self.assertNotAlmostEqual(coordsI[i][j], coordsII[i][j], delta=0.00001)

        # Simulate more to scramble stuff
        sd.get_bd().optimize( 1000 )
        e4 = sd.get_bd().get_scoring_function().evaluate(False)
        print ("Energy %.2f" % e4)
        coordsIV = self._get_beads_coords( sd )

        # Restoration from first simulation through output protobuf file:
        print("*** Initializing positions from ProtoBuf file " + output_pb1)
        f=open(output_pb1, "rb")
        config= IMP.npctransport.Output()
        config.ParseFromString(f.read())
        bch = RMF.BufferConstHandle( config.rmf_conformation )
        fch= RMF.open_rmf_buffer_read_only( bch )
        sd.initialize_positions_from_rmf( fch )
        e5 = sd.get_bd().get_scoring_function().evaluate(False)
        print("*** After initializing positions from RMF file " + output_rmf1)
        print("Energy %.2f" % e5)
        self.assertAlmostEqual(e1, e5, delta=0.00001)
        self.assertNotAlmostEqual(e1, e4, delta=0.00001)
        coordsV = self._get_beads_coords( sd )
        for i in range( len( coordsI ) ):
            for j in range(3):
                self.assertAlmostEqual(coordsI[i][j], coordsV[i][j], delta=0.00001)
                self.assertNotAlmostEqual(coordsI[i][j], coordsIV[i][j], delta=0.00001)

if __name__ == '__main__':
    IMP.test.main()
