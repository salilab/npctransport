from __future__ import print_function
import IMP.npctransport
import IMP.test
import sys
import math
import test_util
#import read_nups


class Tests(IMP.test.TestCase):

    def _create_cfg_file_with_fg_anchors(self, cfg_file):
        '''
        Create a configuration file 'cfg_file' with anchors for fgs

        returns:
        the anchor coordinates 3D tuples by same order as they were added
        to the config file fgs[0] object
        '''
        config = IMP.npctransport.Configuration()
        IMP.npctransport.set_default_configuration(config)
        config.box_is_on.lower=1
        config.box_side.lower=1000
        config.backbone_k.lower=1.0 # kcal/mol/A^2
        config.excluded_volume_k.lower=2 # kcal/mol/A
        config.is_backbone_harmonic=1
        config.backbone_tau_ns.lower=10.0
        config.time_step_factor.lower=1
        config.dump_interval_ns=0.1
        config.statistics_interval_ns=0.1
        config.output_statistics_interval_ns=0.1
        fgs= IMP.npctransport.add_fg_type(config,
                                          type_name="my_fg",
                                          number_of_beads=2, # 16, # 17 is real number for Nup49 ; 15 for Nup57
                                          number=1,
                                          radius=6,
                                          interactions=1,
                                          rest_length_factor = 1.75)
        # dump to file
        f=open(cfg_file, "wb")
        f.write(config.SerializeToString())
        print(config)
        f.close()



    def test_harmonic_spring_score_through_protobuf(self):
        '''
        Test that FG nups can be anchored properly through protobuf file
        '''

        IMP.set_log_level(IMP.SILENT)
        test_util.test_protobuf_installed(self)
        cfg_file = self.get_tmp_file_name("barak_config.pb")
#        assign_file = self.get_tmp_file_name("barak_assign.pb")
        assign_file="output2.pb"
        self._create_cfg_file_with_fg_anchors( cfg_file )
        print("assigning parameter ranges from config")
        num=IMP.npctransport.assign_ranges( cfg_file, assign_file, 0, False, 10 );
        sd= IMP.npctransport.SimulationData(assign_file, False,
                                            "out2.rmf")
#                                            self.get_tmp_file_name("out.rmf"));
        # verify that anchors remain intact during optimization
        if IMP.get_check_level() >= IMP.USAGE_AND_INTERNAL:
            short_init_factor=0.00001
            print("short position initialization in non-fast mode")
        else:
            short_init_factor=0.1
        IMP.npctransport.initialize_positions(sd,[],False,short_init_factor)
        sim_time_ns=5000
        ns_per_chunk=0.01
        sd.get_bd().optimize(100000)
        sd.get_bd().set_current_time(0.0)
        sd.activate_statistics()
        test_util.optimize_in_chunks(sd, sim_time_ns, ns_per_chunk)



if __name__ == '__main__':
    IMP.test.main()
