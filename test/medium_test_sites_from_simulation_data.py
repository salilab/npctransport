from __future__ import print_function
from IMP.npctransport import *
import IMP
import IMP.test
import IMP.core
import sys
import math
#import read_nups
import test_util

fg_R = 15
diffuser_R = 15
fg_coords = [25,25,25]
kap_type = "my_kap"
fg_type = "my_fg"
interaction_k=3
interaction_range=7.5
DEBUG=0

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
        config.excluded_volume_k.lower=7.5
        if(DEBUG):
            config.dump_interval_ns=1.0
#        config.time_step_factor.lower=0.5
#        config.dump_interval_ns=10;
        fgs= IMP.npctransport.add_fg_type(config,
                                          type_name=fg_type,
                                          number_of_beads=1,
                                          number=1,
                                          radius=fg_R,
                                          interactions=6,
                                          rest_length_factor = 1.5)
#        pos=fgs.anchor_coordinates.add()
#        pos.x=fg_coords[0]
#        pos.y=fg_coords[1]
#        pos.z=fg_coords[2]
        kaps= IMP.npctransport.add_float_type(config,
                                              number=1,
                                              radius=diffuser_R,
                                              interactions=6,
                                              type_name=kap_type,
                                              d_factor=1.0
                                              )
        interactionFG_KAP= IMP.npctransport.add_interaction(config,
                                     name0=fg_type,
                                     name1=kap_type,
                                     interaction_k=interaction_k,
                                     interaction_range=interaction_range,
                                                            range_sigma0_deg=60,
                                                            range_sigma1_deg=60)


        # nonspecifics= IMP.npctransport.add_float_type(config,
        #                                               number=1,
        #                                               radius=diffuser_R, #-1,
        #                                               interactions=0)
        # dump to file
        f=open(cfg_file, "wb")
        f.write(config.SerializeToString())
        f.close()

    def find_close_sites(self, sd, p1, p2, distance_thresh=5.0):
        """
        return true if two sites are below distance_thresh from each other
        between particles p1 and p2
        """
        sites1=sd.get_site_centers(IMP.core.Typed(p1).get_type())
        sites2=sd.get_site_centers(IMP.core.Typed(p2).get_type())
        rf1 = IMP.core.RigidBody(p1).get_reference_frame()
        rf2 = IMP.core.RigidBody(p2).get_reference_frame()
        for site1 in sites1:
            site1_c = rf1.get_transformation_to().get_transformed(site1)
            for site2 in sites2:
                site2_c = rf2.get_transformation_to().get_transformed(site2)
                D = IMP.algebra.get_distance(site1_c, site2_c)
                if(D < distance_thresh*1.5 and D >= distance_thresh):
                    print("Almost close sites", site1_c, site2_c, "D", D)
                if(D < distance_thresh):
                    print("Close sites", site1_c, site2_c, "D", D)
                    return True
        return False



    def _assert_kap_in_place(self, sd, really_assert=True):
        first_chain = IMP.npctransport.get_fg_chain(sd.get_fg_chain_roots()[0])
#        print (first_chain)
#        print (vars(first_chain))
        fg_anchor = first_chain.get_bead(0)
#        print (fg_anchor)
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
            self.assertAlmostEqual( D / 15.0, # close up to 15A difference from touching
                                   (fg_R+diffuser_R) / 15.0,
                                   0)
            print("Balls are close", D, " radii: ", fg_R, diffuser_R)
            print ("Asserting close balls")
            self.assertTrue(self.find_close_sites(sd, kap, fg_anchor,
                                                  interaction_range*1.5))
            print ("Asserted close sites")
            print("Kap coords", kap_c, end=' ')
            print("FG coords", anchor_c, end=' ')
            print("D", D)

    def is_stats_interact_(self, output_file):
        ''' verify that stats order parameters know about the interaction '''
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
        in the context of simulation data optimization
        '''

        # Prepare run:
        test_util.test_protobuf_installed(self)
        if IMP.get_check_level() >= IMP.USAGE_AND_INTERNAL:
            print("SLOW MODE")
            fast = False
            short_init_factor=0.0001
            opt_cycles_ns=0.05
            n_iter = 10
            n_good_thresh=1
        else:
            print("FAST MODE")
            fast = True
            short_init_factor=0.01
            opt_cycles_ns=10.0
            n_iter = 100
            n_good_thresh=3
        if(DEBUG):
            opt_cycles_ns=opt_cycles_ns*200
            n_iter=1
            cfg_file = "barak_config.pb"
            assign_file = "barak_assign.pb"
            pymol_file = "sites.pym"
            n_good_thresh=1
        else:
            cfg_file = self.get_tmp_file_name("barak_config.pb")
            assign_file = self.get_tmp_file_name("barak_assign.pb")
            pymol_file = self.get_tmp_file_name("sites.pym")
        self._create_cfg_file_with_fg_anchors( cfg_file )
#        print("assigning parameter ranges from config", cfg_file, end=' ')
#        print("to file", assign_file)
        num=assign_ranges( cfg_file, assign_file, 0, False, 10 );
        sd= IMP.npctransport.SimulationData(assign_file, False)
        if(DEBUG):
            rmf_file="tmp.rmf"
        else:
            rmf_file = self.get_tmp_file_name("out.rmf");
        print("RMF file", rmf_file)
        sd.set_rmf_file(rmf_file, False)
        self._assert_kap_in_place(sd, False)

        # Init and run:
        IMP.set_log_level(IMP.SILENT)
        sd.get_bd().set_log_level(IMP.SILENT)
        IMP.npctransport.initialize_positions( sd, [], False,
                                               short_init_factor)
        opt_cycles_frames=math.ceil(opt_cycles_ns*1E+6/sd.get_bd().get_maximum_time_step());
        print()
        print()
#        IMP.set_log_level(IMP.PROGRESS)
        sd.write_geometry(pymol_file)
        n_good=0
        timer= IMP.npctransport.timer();
        sd.get_bd().set_current_time(0.0)
        sd.activate_statistics()
        sd.get_statistics().reset_statistics_optimizer_states()
        for i in range(n_iter):
            sd.get_bd().optimize(opt_cycles_frames)
            sd.get_statistics().update(timer,opt_cycles_frames)
            try:
                self._assert_kap_in_place(sd, True)
                print ("Asserted kap interacts")
                print ("total energy", sd.get_bd().get_scoring_function().evaluate(False),)
                ppr = sd.get_scoring().get_predicates_pair_restraint()
                print ("predr", ppr.evaluate(False))
                self.assertLess(ppr.evaluate(False), -30.0)
                n_good=n_good+1
                print("NGOOD", n_good)
                self.assertTrue(self.is_stats_interact_(assign_file))
                print("stats interact asserted")
            except AssertionError:
                continue
            if(n_good >= n_good_thresh):
                print("total energy", sd.get_bd().get_scoring_function().evaluate(False), end=' ')
                print("predr", sd.get_scoring().get_predicates_pair_restraint().evaluate(False), end=' ')
                return True
        if fast:
            print("Failed to glue particles after %d iterations x %d opt frames (total %.3f ns)" \
                % (n_iter, opt_cycles_frames, opt_cycles_ns*n_iter))
            self.assertGreaterEqual(n_good, n_good_thresh)
        else:
            print("Debug mode - couldn't glue particles in such short run")


if __name__ == '__main__':
    IMP.test.main()
