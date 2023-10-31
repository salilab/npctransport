from __future__ import print_function
import IMP
import IMP.test
import IMP.npctransport
import IMP.container
import math
import test_util
import glob
import RMF
import numpy as np
import subprocess

radius=8
N_OUTPUT_FRAMES = 3
BOX_SIDE=800

def get_95_conf(rate,time):
    if time==0.0:
        return -1.0
    return 1.96*math.sqrt(rate/time)

class Tests(IMP.test.TestCase):
    
    def _make_cfg(self, cfg_file, is_orientational, is_multiple_hdf5s):
        cfg=test_util.make_simple_cfg(is_slab_on=False, n_particles_factor=1, is_obstacles = 1,
                                      fg_name = "fg", kap_name = "kap", inert_name="inert")
        cfg.box_side.lower=BOX_SIDE
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
        cfg.is_xyz_hist_stats = True
        cfg.is_multiple_hdf5s = is_multiple_hdf5s
        cfg.output_statistics_interval_ns = 1.0
        cfg.simulation_time_ns = N_OUTPUT_FRAMES * cfg.output_statistics_interval_ns
        test_util.write_config_file(cfg_file, cfg)

    def _make_sd(self, is_orientational=False, is_multiple_hdf5s=False):
        cfg_file = self.get_tmp_file_name("my_config.pb")
        assign_file = self.get_tmp_file_name("my_assign.pb")
        self._make_cfg(cfg_file, is_orientational, is_multiple_hdf5s)
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
        n_fgs= G_fgs.get_number_of_children()
        self.assertEqual(n_fgs, 1)
        for j in range(n_fgs):
            fg_type= G_fgs.get_child_name(j)
            self.assertEqual(fg_type, "fg")
            dataset= G_fgs.get_child_int_data_set_3d(fg_type)
            n_voxels_expected = min(500, BOX_SIDE/2)/10 + 5*2 # adding 5 voxels on each side used for slack
            self.assertEqual(dataset.get_size(), 
                             [n_voxels_expected, n_voxels_expected, n_voxels_expected])
            XYZ_as_vector=np.array(dataset.get_block([0,0,0],dataset.get_size()))
        G_floaters= F.get_child_group("floater_xyz_hist")
        n_floaters= G_floaters.get_number_of_children()
        self.assertEqual(n_floaters, 2)
        for j in range(n_floaters):
            floater_type= G_floaters.get_child_name(j)
            print("Floater type: {}".format(floater_type))
            dataset= G_floaters.get_child_int_data_set_3d(floater_type)
            XYZ_as_vector=np.array(dataset.get_block([0,0,0],dataset.get_size()))
                         
    def _test_xyz_stats_using_main(self, 
                                   short_sim_factor,
                                   is_orientational=False):
        """
        test statistics of interaction between an FG and a kap using main.h functions
        """
        for is_multiple_hdf5s in [False, True]:
            print("*** Multiple HDF5s ***" if is_multiple_hdf5s else "*** Single HDF5 ***")
            cfg_file = self.get_tmp_file_name("config{}.pb".format(is_multiple_hdf5s))
            out_file = self.get_tmp_file_name("output{}.pb".format(is_multiple_hdf5s))
            self._make_cfg(cfg_file,
                           is_orientational=is_orientational, 
                           is_multiple_hdf5s=is_multiple_hdf5s)
            # create a temporary python file and call it
            python_file = self.get_tmp_file_name("fake_main.py")
            with open(python_file, 'w') as f:
                f.write(f"""
import IMP.npctransport
cfg_file = "{cfg_file}"
out_file = "{out_file}"
sd = IMP.npctransport.startup(['fake_main', 
    '--configuration', cfg_file, '--output', out_file, 
    '--short_sim_factor', '{short_sim_factor}', '--short_init_factor', '{short_sim_factor}'])
IMP.npctransport.do_main_loop(sd, [])
""")
            subprocess.check_call(["python", python_file])
            hdf5_filenames = glob.glob(out_file + "*.hdf5")
            if is_multiple_hdf5s:
                self.assertEqual(len(hdf5_filenames), np.ceil(N_OUTPUT_FRAMES * short_sim_factor))
            else:
                self.assertEqual(len(hdf5_filenames), 1)
            print("HDF5 files: {}".format(hdf5_filenames))
            for hdf5_filename in hdf5_filenames:
                self._validate_hdf5_file(hdf5_filename)
                    

    def test_hdf5_statistics(self):
        if IMP.get_check_level() >= IMP.USAGE_AND_INTERNAL:
            print("INTERNAL MODE - reduced number of cycles")
            short_sim_factor= 0.01
            print("Short sim factor: {}".format(short_sim_factor))
        else:
            print("FAST MODE")
            short_sim_factor= 1.0
        self._test_xyz_stats_using_main(
            short_sim_factor,
            is_orientational=True)





if __name__ == '__main__':
    IMP.test.main()
