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
import sys


radius=8
N_OUTPUT_FRAMES = 15
N_FULL_OUTPUT_PER_HDF5 = 5
BOX_SIDE_A=300
XYZ_STATS_CROP_FACTOR=1.0
XYZ_STATS_MAX_BOX_SIZE_A=1000
VOXEL_SIZE_A=10
N_FG_BEADS=2
N_FGs=1
N_KAPs=1
N_INERTs=1
N_STATS_PER_HDF5=30.0

def my_make_simple_cfg(outfile, 
                       fg_name = "fg", 
                       kap_name = "kap", 
                       inert_name="inert",
                       is_multiple_hdf5s=False):
    cfg= test_util.get_basic_config()
    IMP.npctransport.add_fg_type(cfg,
                                 type_name=fg_name,
                                 number_of_beads=N_FG_BEADS,
                                 number=N_FGs,
                                 radius=8,
                                 interactions=1,
                                 rest_length_factor = 1.5)
    IMP.npctransport.add_float_type(cfg,
                                    type_name=kap_name,
                                    number=N_KAPs,
                                    radius=20,
                                    interactions=6)
    IMP.npctransport.add_float_type(cfg,
                                    type_name=inert_name,
                                    number=N_INERTs,
                                    radius=20,
                                    interactions=0)
    cfg.box_is_on.lower=1
    cfg.box_side.lower=BOX_SIDE_A
    cfg.slab_is_on.lower= False
    cfg.nonspecific_range.lower= 5.0
    cfg.nonspecific_k.lower= 0.01
    cfg.time_step_factor.lower=8
    cfg.is_xyz_hist_stats=True
    cfg.xyz_stats_crop_factor=XYZ_STATS_CROP_FACTOR
    cfg.xyz_stats_max_box_size_a=XYZ_STATS_MAX_BOX_SIZE_A
    cfg.xyz_stats_voxel_size_a=VOXEL_SIZE_A
    cfg.is_multiple_hdf5s = is_multiple_hdf5s
    cfg.dump_interval_ns=1.0
    cfg.output_statistics_interval_ns = 1.0
    cfg.statistics_interval_ns = cfg.output_statistics_interval_ns / N_STATS_PER_HDF5
    cfg.simulation_time_ns = N_OUTPUT_FRAMES * cfg.output_statistics_interval_ns
    cfg.full_output_statistics_interval_factor = N_FULL_OUTPUT_PER_HDF5
    cfg.statistics_fraction.lower=1.0
    if outfile is not None:
        test_util.write_config_file(outfile, cfg)
    return cfg

class Tests(IMP.test.TestCase):

    def _validate_hdf5_file(self, hdf5_filename, is_multiple_hdf5s):
        print("Handling {}".format(hdf5_filename))
        F = RMF.HDF5.open_file(hdf5_filename)
        G_fgs= F.get_child_group("fg_xyz_hist")
        n_fgs= G_fgs.get_number_of_children()
        self.assertEqual(n_fgs, 1)
        for j in range(n_fgs):
            fg_type= G_fgs.get_child_name(j)
            self.assertEqual(fg_type, "fg")
            dataset= G_fgs.get_child_int_data_set_3d(fg_type)
            n_voxels_expected = XYZ_STATS_CROP_FACTOR * min(XYZ_STATS_MAX_BOX_SIZE_A, BOX_SIDE_A) / VOXEL_SIZE_A \
                                    + 5*2 # adding 5 voxels on each side used for slack
            self.assertEqual(dataset.get_size(), 
                             [n_voxels_expected, n_voxels_expected, n_voxels_expected])
            XYZ_as_vector=np.array(dataset.get_block([0,0,0],dataset.get_size()))
            self.assertEqual(XYZ_as_vector.shape, (n_voxels_expected**3,))
            print("Sum of XYZ: {}".format(np.sum(XYZ_as_vector)))
            sum_observed = np.sum(XYZ_as_vector)
            sum_expected = ( N_FG_BEADS * N_FGs * N_STATS_PER_HDF5  
                              * (1 if is_multiple_hdf5s else N_OUTPUT_FRAMES))
            self.assertAlmostEqual(sum_observed, sum_expected, delta=0.1*sum_observed)
        G_floaters= F.get_child_group("floater_xyz_hist")
        n_floaters= G_floaters.get_number_of_children()
        self.assertEqual(n_floaters, 2)
        for j in range(n_floaters):
            floater_type= G_floaters.get_child_name(j)
            print("Floater type: {}".format(floater_type))
            dataset= G_floaters.get_child_int_data_set_3d(floater_type)
            XYZ_as_vector=np.array(dataset.get_block([0,0,0],dataset.get_size()))
            self.assertEqual(XYZ_as_vector.shape, (n_voxels_expected**3,))
            print("Sum of XYZ: {}".format(np.sum(XYZ_as_vector)))
            sum_observed = np.sum(XYZ_as_vector)
            sum_expected = ( (N_KAPs if floater_type=="kap" else N_INERTs) * N_STATS_PER_HDF5  
                              * (1 if is_multiple_hdf5s else N_OUTPUT_FRAMES))
            self.assertAlmostEqual(sum_observed, sum_expected, delta=0.1*sum_observed)
            
                         
    def _test_xyz_stats_using_main(self, 
                                   short_sim_factor,
                                   is_orientational=False):
        """
        test statistics of interaction between an FG and a kap using main.h functions
        """
        for is_multiple_hdf5s in [False, True]:
            print("*** Multiple HDF5s ***" if is_multiple_hdf5s else "*** Single HDF5 ***")
            cfg_file = "config{}.pb".format(is_multiple_hdf5s)
            out_file = "output{}.pb".format(is_multiple_hdf5s)
            #cfg_file = self.get_tmp_file_name(cfg_file)
            #out_file = self.get_tmp_file_name(out_file)
            my_make_simple_cfg(cfg_file,
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
            subprocess.check_call([sys.executable, python_file])
            hdf5_filenames = glob.glob(out_file + "*.hdf5")
            if is_multiple_hdf5s:
                self.assertEqual(len(hdf5_filenames), np.ceil(N_OUTPUT_FRAMES * short_sim_factor))
            else:
                self.assertEqual(len(hdf5_filenames), 1)
            print("HDF5 files: {}".format(hdf5_filenames))
            for hdf5_filename in hdf5_filenames:
                self._validate_hdf5_file(hdf5_filename, 
                                         is_multiple_hdf5s=is_multiple_hdf5s)
                    

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
