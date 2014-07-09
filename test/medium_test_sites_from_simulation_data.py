from __future__ import print_function
from IMP.npctransport import *
import IMP.base
import IMP.test
import IMP.core
import sys
import math
#import read_nups

fg_R = 25
diffuser_R = 25
fg_coords = [25,25,25]
kap_type = "my_kap"
fg_type = "my_fg"

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
        config.box_side.lower=125
        fgs= IMP.npctransport.add_fg_type(config,
                                          type_name=fg_type,
                                          number_of_beads=1,
                                          number=1,
                                          radius=fg_R,
                                          interactions=12,
                                          rest_length_factor = 1.5)
        pos=fgs.anchor_coordinates.add()
        pos.x=fg_coords[0]
        pos.y=fg_coords[1]
        pos.z=fg_coords[2]
        kaps= IMP.npctransport.add_float_type(config,
                                              number=1,
                                              radius=diffuser_R,
                                              interactions=12,
                                              type_name=kap_type,
                                              d_factor=3
                                              )
        interactionFG_KAP= IMP.npctransport.add_interaction(config,
                                     name0=fg_type,
                                     name1=kap_type,
                                     interaction_k=10,
                                     interaction_range=5)

        # nonspecifics= IMP.npctransport.add_float_type(config,
        #                                               number=1,
        #                                               radius=diffuser_R, #-1,
        #                                               interactions=0)
        # dump to file
        f=open(cfg_file, "wb")
        f.write(config.SerializeToString())
        f.close()

    def find_close_sites(self, sd, p1, p2, distance_thresh=2.0):
        """
        return true if two sites are below distance_thresh from each other
        between particles p1 and p2
        """

        sites1=sd.get_sites(IMP.core.Typed(p1).get_type())
        sites2=sd.get_sites(IMP.core.Typed(p2).get_type())
        rf1 = IMP.core.RigidBody(p1).get_reference_frame()
        rf2 = IMP.core.RigidBody(p2).get_reference_frame()
        for site1 in sites1:
            site1_c = rf1.get_transformation_to().get_transformed(site1)
            for site2 in sites2:
                site2_c = rf2.get_transformation_to().get_transformed(site2)
                D = IMP.algebra.get_distance(site1_c, site2_c)
                if(D < distance_thresh):
                    print("Close sites", site1_c, site2_c, "D", D)
                    return True
        return False



    def _assert_kap_in_place(self, sd, really_assert=True):
        first_chain = IMP.npctransport.get_fg_chain(sd.get_fg_chain_roots()[0])
        print first_chain
        print vars(first_chain)
        print "WHOOO"
        fg_anchor = first_chain.get_bead(0)
        print fg_anchor
        kap = None
        r=sd.get_root()
        for rchild in r.get_children():
            if rchild.get_name()==kap_type:
                kap = rchild.get_child(0)
        assert(kap is not None)
        # check if sphere coordinates are close
        kap_c= IMP.core.XYZ(kap).get_coordinates()
        anchor_c = IMP.core.XYZ( fg_anchor ).get_coordinates()
        D = IMP.algebra.get_distance(kap_c, anchor_c)
        if really_assert:
            # assert that they're close and have close sites
            self.assertAlmostEqual( D / 10.0,
                                   (fg_R+diffuser_R) / 10.0,
                                   0)
            self.assert_( self.find_close_sites(sd, kap, fg_anchor) )
            print("Kap coords", kap_c, end=' ')
            print("FG coords", anchor_c, end=' ')
            print("D", D)

    def is_stats_interact_(self, output_file):
        ''' verify that stats order parametrs know about the interaction '''
        f=open(output_file, "rb")
        o=Output()
        o.ParseFromString(f.read())
        assert(o.statistics is not None)
        for f in o.statistics.floaters:
            if(f.type == kap_type):
                if f.order_params[-1].site_interactions_per_floater > 0.0:
                    return True
        return False



    def test_sites_from_simulation_data(self):
        '''
        Test that the site interaction glues particles together
        in the context of simualtion data optimization
        '''

        if IMP.base.get_check_level() >= IMP.base.USAGE_AND_INTERNAL:
            print("SLOW MODE")
            fast = False
            short_init_factor=0.0001
            opt_cycles=100
            n_iter = 10
            n_good_thresh=1
        else:
            print("FAST MODE")
            fast = True
            short_init_factor=0.1
            opt_cycles=12500
            n_iter = 200
            n_good_thresh=3

        # prepare run
        cfg_file = self.get_tmp_file_name("barak_config.pb")
        assign_file = self.get_tmp_file_name("barak_assign.pb")
        pymol_file = self.get_tmp_file_name("sites.pym")
        self._create_cfg_file_with_fg_anchors( cfg_file )
        print("assigning parameter ranges from config", cfg_file, end=' ')
        print("to file", assign_file)
        num=assign_ranges( cfg_file, assign_file, 0, False, 10 );
        rmf_file = self.get_tmp_file_name("out.rmf");
        print("RMF file", rmf_file)
        sd= IMP.npctransport.SimulationData(assign_file, False)
        sd.set_rmf_file(rmf_file, False)
        self._assert_kap_in_place(sd, False)
        # init and run
        IMP.base.set_log_level(IMP.base.SILENT)
        sd.get_bd().set_log_level(IMP.base.SILENT)
        IMP.npctransport.initialize_positions( sd, [], False,
                                               short_init_factor)
        print()
        print()
#        IMP.base.set_log_level(IMP.base.PROGRESS)
        sd.write_geometry(pymol_file)
        n_good=0
        timer= IMP.npctransport.timer();
        sd.get_statistics().reset_statistics_optimizer_states()
        sd.get_bd().set_current_time(0.0)
        for i in range(n_iter):
            sd.get_bd().optimize(opt_cycles)
            sd.get_statistics().update(timer,opt_cycles)
            try:
                self._assert_kap_in_place(sd, True)
                print "total energy", sd.get_bd().get_scoring_function().evaluate(False),
                print "predr", sd.get_scoring().get_predicates_pair_restraint().evaluate(False)
                self.assert_(sd.get_scoring().get_predicates_pair_restraint().evaluate(False) < -30.0)
                self.assert_(self.is_stats_interact_(assign_file))
                n_good=n_good+1
            except:
                continue
            if(n_good >= n_good_thresh):
                print("total energy", sd.get_bd().get_scoring_function().evaluate(False), end=' ')
                print("predr", sd.get_scoring().get_predicates_pair_restraint().evaluate(False), end=' ')
                return True
        if fast:
            print("Failed to glue particles after %d iterations of %d opt cycles" \
                % (n_iter, opt_cycles))
            self.assert_(n_good >= n_good_thresh)
        else:
            print("Debug mode - couldn't glue particles in such short run")


if __name__ == '__main__':
    IMP.test.main()
