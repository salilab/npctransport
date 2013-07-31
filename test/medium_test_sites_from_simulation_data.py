#!/usr/bin/python
from IMP.npctransport import *
import IMP.base
import IMP.test
import sys
import math
#import read_nups

fg_R = 25
diffuser_R = 25
fg_coords = [25,25,25]

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
        config.box_is_on.lower=1
        config.box_side.lower=150
        fgs= IMP.npctransport.add_fg_type(config,
                                          type_name="my_fg",
                                          number_of_beads=1,
                                          number=1,
                                          radius=fg_R,
                                          interactions=7,
                                          rest_length_factor = 1.5)
        coords = []
        pos=fgs.anchor_coordinates.add()
        c = fg_coords
        pos.x=c[0]
        pos.y=c[1]
        pos.z=c[2]
        kaps= IMP.npctransport.add_float_type(config,
                                              number=1,
                                              radius=diffuser_R,
                                              interactions=12,
                                              type_name="my_kap",
                                              d_factor=3
                                              )
        interactionFG_KAP= IMP.npctransport.add_interaction(config,
                                     name0="my_fg",
                                     name1="my_kap",
                                     interaction_k=10,
                                     interaction_range=10)

        # nonspecifics= IMP.npctransport.add_float_type(config,
        #                                               number=1,
        #                                               radius=diffuser_R, #-1,
        #                                               interactions=0)
        # dump to file
        f=open(cfg_file, "wb")
        f.write(config.SerializeToString())
        print config
        f.close()

    def _assert_anchors_in_place(self, sd):
        '''
        verify that fg anchors are still in same place as in the originally
        specified coordinates coords
        '''
        fgs = sd.get_fg_chains( ) # atom.Hierarchies
#        print " fgs ", fgs
        # verify that anchors are in placed
        self.assert_(len(fgs)==1)
        fg=fgs[0]
        anchor_c = IMP.core.XYZ( fg.get_child(0) ).get_coordinates()
#        print "fg_sites", sd.get_sites(IMP.core.Typed(fg).get_type())
        d2 = sum([(a-b)**2 for a,b in zip(anchor_c,fg_coords)])
#        print "cur fg anchor / orig anchor / dist2", anchor_c, fg_coords, d2
        self.assertEqual(d2, 0)

    def _assert_kap_in_place(self, sd, really_assert=True):
        r=sd.get_root()
        for p in r.get_children():
            if p.get_name()=="my_kap":
                for kap in p.get_children():
                    kap_c= IMP.core.XYZ(kap).get_coordinates()
                    print "kap ", kap_c,
#                    print "kap sites",
#                    print sd.get_sites(IMP.core.Typed(kap).get_type())
                    d2 = sum([(x-y)**2 for x,y in zip(kap_c,fg_coords)])
#                    print d2
                    print "d=", math.sqrt(d2)
                    if really_assert:
                        self.assertAlmostEqual(d2/100.0,
                                               (fg_R+diffuser_R)**2 / 100,
                                               0)

    def test_sites_from_simulation_data(self):
        '''
        Test that the site interaction glues particles together
        in the context of simualtion data optimization
        '''

        IMP.base.set_log_level(IMP.SILENT)
        cfg_file = self.get_tmp_file_name("barak_config.pb")
        assign_file = self.get_tmp_file_name("barak_assign.pb")
        self._create_cfg_file_with_fg_anchors( cfg_file )
        print "assigning parameter ranges from config"
        num=assign_ranges( cfg_file, assign_file, 0, False, 10 );
        sd= IMP.npctransport.SimulationData(assign_file, False,
                                            self.get_tmp_file_name("out.rmf"));
        self._assert_anchors_in_place(sd)
        # verify that anchors remain intact during optimization
#        if IMP.build == "fast" or IMP.build == "release":
        print "Check level ", IMP.base.get_check_level()
        print "Usage and internal: ", IMP.base.USAGE_AND_INTERNAL
        if IMP.base.get_check_level() >= IMP.base.USAGE_AND_INTERNAL:
            print "SLOW MODE"
            fast = False
            short_init_factor=0.0001
            opt_cycles=100
            n_iter = 10
        else:
            print "FAST MODE"
            fast = True
            short_init_factor=0.1
            opt_cycles=25000
            n_iter = 100
        self._assert_kap_in_place(sd, False)
        IMP.npctransport.initialize_positions( sd, [], False,
                                               short_init_factor)
        self._assert_anchors_in_place(sd)
        self._assert_kap_in_place(sd, False)
        print "total", sd.get_bd().get_scoring_function().evaluate(False),
        print "predr", sd.get_scoring().get_predr().evaluate(False)
        for i in range(n_iter):
            sd.get_bd().optimize(opt_cycles)
            print "total", sd.get_bd().get_scoring_function().evaluate(False),
            print "predr", sd.get_scoring().get_predr().evaluate(False),
            try:
                self._assert_anchors_in_place(sd)
                self._assert_kap_in_place(sd, True)
            except:
                continue
            return True
        if fast:
            self._assert_kap_in_place(sd, True)
            self.assert_(sd.get_scoring().get_predr().evaluate(False) < -80.0)
        else:
            print "Debug mode - couldn't glue particles in such short run"


if __name__ == '__main__':
    IMP.test.main()
