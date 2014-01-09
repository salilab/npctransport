#!/usr/bin/python
from IMP.npctransport import *
import IMP.base
import IMP.test
import sys
import math
#import read_nups


class Tests(IMP.test.TestCase):

    def _create_cfg_file_with_fgs(self, cfg_file):
        '''
        Create a configuration file 'cfg_file' with fg chains
        '''
        config = Configuration()
        IMP.npctransport.set_default_configuration(config)
        config.box_is_on.lower=1
        config.box_side.lower=200
        fgs= IMP.npctransport.add_fg_type(config,
                                          type_name="my_fg1",
                                          number_of_beads=3,
                                          number=2,
                                          radius=25,
                                          interactions=7,
                                          rest_length_factor = 1.5)
        fgs= IMP.npctransport.add_fg_type(config,
                                          type_name="my_fg2",
                                          number_of_beads=6,
                                          number=3,
                                          radius=25,
                                          interactions=7,
                                          rest_length_factor = 1.5)
        # diffuser_R=25
        # kaps= IMP.npctransport.add_float_type(config,
        #                                       number=1,
        #                                       radius=diffuser_R,
        #                                       interactions=12,
        #                                       type_name="my_kap",
        #                                       d_factor=3
        #                                       )

        # nonspecifics= IMP.npctransport.add_float_type(config,
        #                                               number=1,
        #                                               radius=diffuser_R, #-1,
        #                                               interactions=0)
        # dump to file
        f=open(cfg_file, "wb")
        f.write(config.SerializeToString())
#        print config
        f.close()


    def test_get_fg_chains(self):
        '''
        Test that SimulationData.get_fg_chains() works correctly
        '''

        IMP.base.set_log_level(IMP.SILENT)
        cfg_file = self.get_tmp_file_name("barak_config.pb")
        assign_file = self.get_tmp_file_name("barak_assign.pb")
        coords = self._create_cfg_file_with_fgs(cfg_file)
        print "assigning parameter ranges from config"
        num=assign_ranges( cfg_file, assign_file, 0, False, 10 )
        sd= IMP.npctransport.SimulationData(assign_file, False)
        sd.set_rmf_file( self.get_tmp_file_name("out.rmf") )
        fgs = sd.get_fg_chains( ) # atom.Hierarchies
        print " fgs ", fgs
        ones = 0
        twos = 0
        for fg in fgs:
            type_name = IMP.core.Typed(fg.get_child(0)).get_type().get_string()
            if(type_name == "my_fg1"):
                self.assert_(fg.get_number_of_children() == 3)
                ones = ones + 1
            if(type_name == "my_fg2"):
                self.assert_(fg.get_number_of_children() == 6)
                twos = twos + 1
        self.assert_(ones == 2 and twos == 3)


if __name__ == '__main__':
    IMP.test.main()
