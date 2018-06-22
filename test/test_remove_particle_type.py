from __future__ import print_function
import IMP
import IMP.test
import IMP.npctransport
import IMP.container
import math
import test_util
import sys

radius=8

def get_95_conf(rate,time):
    if time==0.0:
        return -1.0
    return 1.96*math.sqrt(rate/time)

class Tests(IMP.test.TestCase):

    def _make_sd(self, is_orientational=False, extra_fg_types=['my_fg2']):
        cfg_file = self.get_tmp_file_name("barak_config.pb")
        assign_file = self.get_tmp_file_name("barak_assign.pb")
        cfg=test_util.make_simple_cfg(is_slab_on=False, n_particles_factor=1)
        for extra_fg_type in extra_fg_types:
            fg2= IMP.npctransport.add_fg_type(cfg,
                                              type_name=extra_fg_type,
                                              number_of_beads=2,
                                              number=1,
                                              radius=6,
                                              interactions=1,
                                              rest_length_factor = 1.5)
        cfg.box_side.lower=350
        cfg.time_step_factor.lower=8
        cfg.statistics_interval_ns=0.01
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



    def _get_particles_of_type(self, sd, pt):
        """
        get bead particles of type pt in simulation data sd
        """
        ps=sd.get_beads()
        ret=[]
        for cur_p in ps:
            cur_pt=IMP.core.Typed(cur_p).get_type()
            if(cur_pt==pt):
                ret.append(cur_p)
        return ret


    def test_remove_particle_type(self):
        pt1=IMP.core.ParticleType("my_fg")
        pt2=IMP.core.ParticleType("my_fg2")
        sd= self._make_sd(False, ["my_fg2"])
        timer = IMP.npctransport.create_boost_timer()
        sd.activate_statistics()
#        sd.set_rmf_file("1.rmf")
        self._run_sd(sd, 50)
        self.assertEqual(len(self._get_particles_of_type(sd, pt1)), 2)
        self.assertEqual(len(self._get_particles_of_type(sd, pt2)), 2)
        sd.remove_particle_type(pt2)
#        sd.set_rmf_file("2.rmf")
        self._run_sd(sd, 50)
        self.assertEqual(len(self._get_particles_of_type(sd, pt1)), 2)
        self.assertEqual(len(self._get_particles_of_type(sd, pt2)), 0)

    def test_remove_particle_typ(self):
        extra_fg_types_strings=["Nup1", "Nup1b", "Nup1c", "Nup159", "Nup2"]
        fg_types= [IMP.core.ParticleType(s_fg_type) for s_fg_type in extra_fg_types_strings]
        fg_types.append(IMP.core.ParticleType("my_fg"))
        sd= self._make_sd(False, extra_fg_types_strings)
        timer = IMP.npctransport.create_boost_timer()
        sd.activate_statistics()
#        sd.set_rmf_file("1.rmf")
        self._run_sd(sd, 50)
        for fg_type in fg_types:
            self.assertEqual(len(self._get_particles_of_type(sd, fg_type)), 2)
        sd.remove_fgs_with_prefix("Nup1")
#        sd.set_rmf_file("2.rmf")
        self._run_sd(sd, 50)
        correct_removed_indexes=[0,1,2]
        correct_kept_indexes=[3,4,5]
        for i in correct_removed_indexes:
            self.assertEqual(len(self._get_particles_of_type(sd, fg_types[i])),
                             0)
        for i in correct_kept_indexes:
            self.assertEqual(len(self._get_particles_of_type(sd, fg_types[i])),
                             2)


if __name__ == '__main__':
    IMP.test.main()
