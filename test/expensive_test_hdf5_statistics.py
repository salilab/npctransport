from __future__ import print_function
import IMP
import IMP.test
import IMP.npctransport
import IMP.container
import math
import test_util
import glob
import RMF

radius=8

def get_95_conf(rate,time):
    if time==0.0:
        return -1.0
    return 1.96*math.sqrt(rate/time)

class Tests(IMP.test.TestCase):

    def _make_sd(self, is_orientational=False, is_multiple_hdf5s=False):
        cfg_file = self.get_tmp_file_name("my_config.pb")
        assign_file = self.get_tmp_file_name("my_assign.pb")
        print("Configuration file: " + cfg_file)
        print("Assignment file: " + assign_file)
        cfg=test_util.make_simple_cfg(is_slab_on=False, n_particles_factor=1, is_obstacles = 1,
                                      fg_name = "fg", kap_name = "kap", inert_name="inert")
        cfg.box_side.lower=35
        cfg.time_step_factor.lower=8
        cfg.statistics_interval_ns=0.01
        for i in cfg.interactions:
            if i.type0=="kap0" or i.type1=="kap":
                i.interaction_range.lower=10
                if(is_orientational):
                    i.range_sigma0_deg.lower=45
                    i.range_sigma1_deg.lower=45
                    i.interaction_k.lower=0.18
                else:
                    i.interaction_k.lower=0.45
        cfg.is_multiple_hdf5s = is_multiple_hdf5s
        cfg.output_statistics_interval_ns = 50.0
        test_util.write_config_file(cfg_file, cfg)
        num=IMP.npctransport.assign_ranges( cfg_file, assign_file, 0,
                           False, 10 );
        sd= IMP.npctransport.SimulationData(assign_file, False)
        return sd

    def _run_sd(self, sd, n_cycles):
        self.assertTrue(sd != None)
        IMP.set_log_level(IMP.SILENT)
        sd.get_bd().set_log_level(IMP.SILENT)
        time_step_fs=sd.get_bd().get_maximum_time_step()
        print("Simulating for", n_cycles*time_step_fs*1e-6, "ns")
        sd.get_bd().optimize(n_cycles)
        print(sd.get_bd().get_scoring_function().evaluate(False))
        sd.get_statistics().update(IMP.npctransport.create_boost_timer(),
                                   n_cycles)

    def _validate_hdf5_file(self, hdf5_filename):
        print("Handling {}".format(hdf5_filename))
        F = RMF.HDF5.open_file(hdf5_filename)
        G_fgs= F.get_child_group("fg_xyz_hist")
        G_floaters= F.get_child_group("floater_xyz_hist")
        

    def _process_xyz_stat(self, sd, is_orientational=False):
        output_filename=sd.get_statistics().get_output_file_name()
        hdf5_filenames = glob.glob(output_filename + "*.hdf5")
        for hdf5_filename in hdf5_filenames:
            self._validate_hdf5_file(hdf5_filename)

    def _test_xyz_stats(self, n_cycles, short_init_factor,
                                n_trials, is_orientational=False):
        """
        test statistics of interaction between an FG and a kap
        """
        for is_multiple_hdf5s in [False, True]:
            print("Multiple HDF5s" if is_multiple_hdf5s else "Single HDF5")
            sd=self._make_sd(is_orientational=is_orientational, 
                             is_multiple_hdf5s=is_multiple_hdf5s)
            IMP.npctransport.initialize_positions( sd, [], False,
                                                short_init_factor)
            sd.get_bd().optimize(n_cycles) # still without stats
            sd.get_bd().set_current_time(0.0);
            sd.get_statistics().reset_statistics_optimizer_states();
            sd.activate_statistics()
            while n_trials>0:
                try:
                    self._run_sd(sd, n_cycles)
                    self._process_xyz_stat(sd, is_orientational)
                    break
                except AssertionError as e:
                    print("Failed -", n_trials, "left")
                    n_trials=n_trials-1
                    if(n_trials==0):
                        raise e

    def test_hdf5_statistics(self):
        if IMP.get_check_level() >= IMP.USAGE_AND_INTERNAL:
            print("INTERNAL MODE - reduced number of cycles")
            n_cycles= 10
            n_trials= 1
            short_init_factor= 0.01
        else:
            print("FAST MODE")
            n_cycles=100000
            n_trials= 200
            short_init_factor= 0.1
        self._test_xyz_stats(n_cycles,
                                     short_init_factor,
                                     n_trials,is_orientational=True)





if __name__ == '__main__':
    IMP.test.main()
