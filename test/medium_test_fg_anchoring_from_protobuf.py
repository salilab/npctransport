from IMP.npctransport import *
import IMP.test
import sys
import math
#import read_nups


class Tests(IMP.test.TestCase):

    def _create_cfg_file_with_fg_anchors(self, cfg_file):
        '''
        Create a configuration file 'cfg_file' with anchors for fgs

        returns:
        the anchor coordinates 3D tuples by same order as they were added
        to the config file fgs[0] object
        '''
        config = Configuration()
        IMP.npctransport.set_default_configuration(config)
        kaps_R = 10.0
        config.box_is_on.lower=1
        config.box_side.lower=300
        fgs= IMP.npctransport.add_fg_type(config,
                                          type_name="my_fg",
                                          number_of_beads=2, # 16, # 17 is real number for Nup49 ; 15 for Nup57
                                          number=3,
                                          radius=6,
                                          interactions=1,
                                          rest_length_factor = 1.5)
        coords = []
        for i in range(10):
            pos=fgs.anchor_coordinates.add()
            c = (i, i*10, i*5)
            pos.x=c[0]
            pos.y=c[1]
            pos.z=c[2]
            coords.append(c)
#        kaps= IMP.npctransport.add_float_type(config,
#                                              number=1,
#                                              radius=kaps_R,
#                                              interactions=12)
#        nonspecifics= IMP.npctransport.add_float_type(config,
#                                                      number=1,
#                                                      radius=kaps_R, #-1,
#                                                      interactions=0)
        # dump to file
        f=open(cfg_file, "wb")
        f.write(config.SerializeToString())
        print config
        f.close()
        return coords

    def _assert_anchors_in_place(self, sd, coords):
        '''
        verify that fg anchors are still in same place as in the originally
        specified coordinates coords
        '''
        fgs = sd.get_fg_chains(  ) # atom.Hierarchies
        # verify that anchors are in placed
        for fg,c in zip(fgs,coords):
            fg_c = IMP.core.XYZ( fg.get_child(0) ).get_coordinates()
            d2 = sum([(a-b)**2 for a,b in zip(fg_c,c)])
            print fg_c, c, d2
            self.assertEqual(d2, 0)

    def test_fg_anchoring_through_protobuf(self):
        '''
        Test that FG nups can be anchored properly through protobuf file
        '''

        IMP.base.set_log_level(IMP.SILENT)
        cfg_file = self.get_tmp_file_name("barak_config.pb")
        assign_file = self.get_tmp_file_name("barak_assign.pb")
        coords = self._create_cfg_file_with_fg_anchors( cfg_file )
        print "assigning parameter ranges from config"
        num=assign_ranges( cfg_file, assign_file, 0, False, 10 );
        sd= IMP.npctransport.SimulationData(assign_file, False, "")
#                                            self.get_tmp_file_name("out.rmf"));
        self._assert_anchors_in_place(sd, coords)
        # verify that anchors remain intact during optimization
        if IMP.base.get_check_level() >= IMP.base.USAGE_AND_INTERNAL:
            short_init_factor=0.00001
            print "short position initialization in non-fast mode"
        else:
            short_init_factor=0.1
        IMP.npctransport.initialize_positions(sd,[],False,short_init_factor)
        self._assert_anchors_in_place(sd, coords)
        sd.get_bd().optimize(1000)
        self._assert_anchors_in_place(sd, coords)



if __name__ == '__main__':
    IMP.test.main()
