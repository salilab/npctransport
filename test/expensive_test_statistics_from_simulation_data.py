from __future__ import print_function
import IMP
import IMP.test
import IMP.npctransport
import IMP.container
import math
import test_util

radius=8

def get_95_conf(rate,time):
    if time==0.0:
        return -1.0
    return 1.96*math.sqrt(rate/time)

class Tests(IMP.test.TestCase):

    def _make_sd(self, is_orientational=False):
        cfg_file = self.get_tmp_file_name("barak_config.pb")
        assign_file = self.get_tmp_file_name("barak_assign.pb")
        assign_file='tmp.pb'
        cfg=test_util.make_simple_cfg(is_slab_on=False, n_particles_factor=1)
        cfg.box_side.lower=35
        cfg.time_step_factor.lower=8
        cfg.statistics_interval_ns=0.01
        for i in cfg.interactions:
            if i.type0=="kap0" or i.type1=="kap0":
                i.interaction_range.lower=10
                if(is_orientational):
                    i.range_sigma0_deg.lower=45
                    i.range_sigma1_deg.lower=45
                    i.interaction_k.lower=0.18
                else:
                    i.interaction_k.lower=0.45
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


    def _process_sd_stat(self, sd, is_orientational=False):
        assign_file=sd.get_statistics().get_output_file_name()
        o= IMP.npctransport.Output()
        o.ParseFromString(open(assign_file,'r').read())
        time_ns=o.statistics.bd_simulation_time_ns
        for i in o.statistics.interactions:
            print(i.type0,i.type1)
            if(i.type0<>"kap0" and i.type1<>"kap0"):
                continue
            n=len(i.order_params)+0.0
            koff_i=0.0
            off_i_time_ns=0.0
            kon_i=0.0
            on_i_time_ns=0.0
            koff_ii=0.0
            off_ii_time_ns=0.0
            kon_ii=0.0
            on_ii_time_ns=0.0
            fb_i=0.0
            fb_ii=0.0
            misc_time_ns=0.0
            for iop in i.order_params:
                # Get params and compute confidence intervals
                koff_i= koff_i + iop.avg_off_per_bound_i_per_ns*iop.off_i_stats_period_ns
                off_i_time_ns= off_i_time_ns + iop.off_i_stats_period_ns
                kon_i= kon_i + iop.avg_on_per_unbound_i_per_ns*iop.on_i_stats_period_ns
                on_i_time_ns= on_i_time_ns + iop.on_i_stats_period_ns
                koff_ii= koff_ii + iop.avg_off_per_bound_ii_per_ns*iop.off_ii_stats_period_ns
                off_ii_time_ns= off_ii_time_ns + iop.off_ii_stats_period_ns
                kon_ii= kon_ii + iop.avg_on_per_unbound_ii_per_ns*iop.on_ii_stats_period_ns
                on_ii_time_ns= on_ii_time_ns + iop.on_ii_stats_period_ns
                fb_i= fb_i + iop.avg_fraction_bound_particles_i*iop.misc_stats_period_ns
                fb_ii= fb_ii + iop.avg_fraction_bound_particles_ii*iop.misc_stats_period_ns
                misc_time_ns= misc_time_ns + iop.misc_stats_period_ns
            koff_i= koff_i/(off_i_time_ns+0.0000001)
            kon_i= kon_i/(on_i_time_ns+0.000001)
            koff_ii= koff_ii/(off_ii_time_ns+0.00001)
            kon_ii= kon_ii/(on_ii_time_ns+0.000001)
            fb_i= fb_i/(misc_time_ns+0.000001)
            fb_ii= fb_ii/(misc_time_ns+0.000001)
            # Get confience intervals
            conf95_koff_i=get_95_conf(koff_i,off_i_time_ns)
            conf95_kon_i=get_95_conf(kon_i,on_i_time_ns)
            conf95_koff_ii=get_95_conf(koff_ii,off_ii_time_ns)
            conf95_kon_ii=get_95_conf(kon_ii,on_ii_time_ns)
            # Print stats
            print("Off_i", koff_i, "+-", conf95_koff_i)
            print("On_i", kon_i, "+-", conf95_kon_i)
            print("Off_ii", koff_ii, "+-", conf95_koff_ii)
            print("On_ii", kon_ii, "+-", conf95_kon_ii)
            print("Kd_i",  koff_i /(kon_i+0.00000001) )
            print("Kd_ii", koff_ii/(kon_ii+0.00000001))
            print("%% bound: I %.1f%% II %.1f%%" % ( 100*fb_i, 100*fb_ii))
 #            # Verify results
            if IMP.get_check_level() >= IMP.USAGE_AND_INTERNAL:
                return
            if is_orientational:
                DELTA_FACTOR=0.5
                self.assertAlmostEqual(koff_i,2.21,delta=0.45*DELTA_FACTOR)
                self.assertAlmostEqual(kon_i,1.15,delta=0.25*DELTA_FACTOR)
                self.assertAlmostEqual(koff_ii,1.76,delta=0.35*DELTA_FACTOR)
                self.assertAlmostEqual(kon_ii,2.16,delta=0.45*DELTA_FACTOR)
                self.assertAlmostEqual(fb_i,0.34,delta=0.11*DELTA_FACTOR)
                self.assertAlmostEqual(fb_ii,0.55,delta=0.15*DELTA_FACTOR)
            else:
                DELTA_FACTOR=0.5
                self.assertAlmostEqual(koff_i,1.81,delta=0.35*DELTA_FACTOR)
                self.assertAlmostEqual(kon_i,0.75,delta=0.2*DELTA_FACTOR)
                self.assertAlmostEqual(koff_ii,1.41,delta=0.3*DELTA_FACTOR)
                self.assertAlmostEqual(kon_ii,1.61,delta=0.33*DELTA_FACTOR)
                self.assertAlmostEqual(fb_i,0.285,delta=0.1*DELTA_FACTOR)
                self.assertAlmostEqual(fb_ii,0.524,delta=0.15*DELTA_FACTOR)
        return;


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

    def _test_interaction_stats(self, n_cycles, short_init_factor,
                                n_trials, is_orientational=False):
        """
        test statistics of interaction between an FG and a kap
        """
        # Non-orientational
        sd=self._make_sd(is_orientational)
        IMP.npctransport.initialize_positions( sd, [], False,
                                               short_init_factor)
        sd.get_bd().optimize(n_cycles) # still without stats
        sd.get_bd().set_current_time(0.0);
        sd.get_statistics().reset_statistics_optimizer_states();
        sd.activate_statistics()
        while n_trials>0:
            try:
                self._run_sd(sd, n_cycles)
                self._process_sd_stat(sd, is_orientational)
                break
            except AssertionError as e:
                print("Failed -", n_trials, "left")
                n_trials=n_trials-1
                if(n_trials==0):
                    raise e

    def test_all_interaction_stats(self):
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
        print("-- Non-orientational --")
        self._test_interaction_stats(n_cycles,short_init_factor,
                                     n_trials,is_orientational=False)
        print("======================")
        print("-- Orientational --")
        self._test_interaction_stats(n_cycles,short_init_factor,
                                     n_trials,is_orientational=True)

    def test_statistics_activation(self):
        sd= self._make_sd(False)
        timer = IMP.npctransport.create_boost_timer()
        self.assertFalse( sd.get_statistics().get_is_activated() )
        with self.assertRaises(IMP.UsageException):
            sd.get_statistics().update( timer )
        sd.activate_statistics()
        self.assertTrue( sd.get_statistics().get_is_activated() )
        sd.get_statistics().update( timer ) # should not throw an exception now




if __name__ == '__main__':
    IMP.test.main()
