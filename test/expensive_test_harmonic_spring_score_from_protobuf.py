from __future__ import print_function
import IMP.npctransport
import IMP.core
import IMP.algebra
import IMP.test
import sys
import math
import test_util

RADIUS= 7
REST_LENGTH_FACTOR= 1.75
TAU_NS= 1.0
NS_PER_CHUNK= 0.5*TAU_NS
BACKBONE_K= 0.5

class Tests(IMP.test.TestCase):

    def _do_particles_report(self, m, pis):
        global RADIUS
        global REST_LENGTH_FACTOR
        for pi in pis:
            print("Particle index",pi)
            if IMP.core.XYZR.get_is_setup(m, pi):
                xyzr= IMP.core.XYZR(m, pi)
                print("XYZR", xyzr)
            if IMP.atom.Diffusion.get_is_setup(m, pi):
                d= IMP.atom.Diffusion(m, pi)
                print("Diffusion", d)
            if IMP.npctransport.RelaxingSpring.get_is_setup(m, pi):
                rs= IMP.npctransport.RelaxingSpring(m, pi)
                print("RelaxingSpring:", rs)
                eq_rest_length_factor= rs.get_equilibrium_rest_length_factor()
                self.assertAlmostEqual(eq_rest_length_factor,
                                       REST_LENGTH_FACTOR, delta=1e-6)
                expected_rest_length= rs.get_equilibrium_rest_length_factor() * IMP.core.XYZR(m,pi).get_radius() * 2.0
                actual_rest_length= rs.get_rest_length()
                print("Equilibrium rest length expected", expected_rest_length, "A")
                print("Actual rest length", actual_rest_length, "A")
                allowed_delta= math.sqrt(3.0/BACKBONE_K) # delta G ~ BACKBONE_K*delta^2 kcal/mol, so bound at 3.0 kcal/mol
                print("Allowed deviation within order of 3 kcal/mol", allowed_delta, "A")
                self.assertAlmostEqual(expected_rest_length,
                                       actual_rest_length,
                                       delta= allowed_delta)



    def _create_cfg_file_with_fg_anchors(self, cfg_file, n_beads=2):
        '''
        Create a configuration file 'cfg_file' with anchors for fgs

        returns:
        the anchor coordinates 3D tuples by same order as they were added
        to the config file fgs[0] object
        '''
        global RADIUS
        global REST_LENGTH_FACTOR
        global TAU_NS
        config = IMP.npctransport.Configuration()
        IMP.npctransport.set_default_configuration(config)
        config.angular_D_factor.lower=1
        config.nonspecific_k.lower=.000075
        config.box_is_on.lower=1
        config.box_side.lower=10000
        config.backbone_k.lower=BACKBONE_K # kcal/mol/A^2
        config.excluded_volume_k.lower=2 # kcal/mol/A
        config.is_backbone_harmonic=1
        config.backbone_tau_ns.lower= TAU_NS
        config.time_step_factor.lower=2
        config.dump_interval_ns=0.1
        config.statistics_interval_ns=0.1
        config.output_statistics_interval_ns=0.1
        config.temperature_k.lower=298
        fgs= IMP.npctransport.add_fg_type(config,
                                          type_name="my_fg",
                                          number_of_beads=n_beads, # 16, # 17 is real number for Nup49 ; 15 for Nup57
                                          number=1,
                                          radius=RADIUS,
                                          interactions=1,
                                          rest_length_factor = REST_LENGTH_FACTOR)
        # dump to file
        f=open(cfg_file, "wb")
        f.write(config.SerializeToString())
        #        print(config)
        f.close()

    def _analyze_run(self, output_file, R):
        ''' R is the rest length distance data from the RMF file '''
        global TAU_NS
        global NS_PER_CHUNK
        f=open(output_file, "rb")
        output= IMP.npctransport.Output()
        output.ParseFromString(f.read())
        stats= output.statistics
        B=[]
        E=[]
        T=[0.0]
        for fg in stats.fgs:
            for op in fg.order_params:
                B.append(op.mean_bond_distance)
                E.append(op.mean_end_to_end_distance)
                T.append(op.time_ns)
                # end-to-end and bond distances ought to be equal for two bead chains
                if output.assignment.fgs[0].number_of_beads.value==2:
                    self.assertAlmostEqual(op.mean_end_to_end_distance,
                                           op.mean_bond_distance, delta=1e-6)
                    self.assertAlmostEqual(op.mean_square_end_to_end_distance,
                                           op.mean_square_bond_distance,
                                           delta=1e-6)
        if not test_util.check_import_pandas_with_series_autocorr():
            print("WARNING: pandas module not installed or too old, skipping autocorrelation test")
            return
        import pandas as pd
        i_tau= int(max(round(TAU_NS/NS_PER_CHUNK), 1))
        max_i= min(10*i_tau, len(E))
        print("time[ns]  corr-particles  corr-bond-rest-length")
        for i in range(max_i):
            mean_ACR=0.0
            for j in range(R.shape[1]):
                mean_ACR= mean_ACR+pd.Series(R[:,j]).autocorr(i) / R.shape[1]
            print("%8.3f  " % T[i], \
                      "%10.3f  " % pd.Series(B).autocorr(i), \
                      "%10.3f  " % pd.Series(E).autocorr(i), \
                      "%10.3f  " % mean_ACR)
        if len(E)>i_tau:
            self.assertAlmostEqual(
                pd.Series(R[:,0]).autocorr(i_tau),
                math.exp(-1),
                delta=0.07)
        if len(E)>2*i_tau:
            self.assertAlmostEqual(
                pd.Series(R[:,0]).autocorr(2*i_tau),
                math.exp(-2),
                delta=0.07)
        if len(E)>5*i_tau:
            self.assertAlmostEqual(
                pd.Series(R[:,0]).autocorr(5*i_tau),
                math.exp(-5),
                delta=0.07)



    def _test_harmonic_spring_score_from_protobuf(self, n_beads=2):
        '''
        Test that FG nups can be anchored properly through protobuf file
        '''
        global NS_PER_CHUNK

        IMP.set_log_level(IMP.SILENT)
        test_util.test_protobuf_installed(self)
        cfg_file = self.get_tmp_file_name("barak_config.pb")
#        assign_file = self.get_tmp_file_name("barak_assign.pb")
        assign_file="output2.pb"
        self._create_cfg_file_with_fg_anchors( cfg_file, n_beads )
        print("assigning parameter ranges from config")
        num=IMP.npctransport.assign_ranges( cfg_file, assign_file, 0, False, 10 );
        sd= IMP.npctransport.SimulationData(assign_file,
                                            False,
                                            "")
#                                            "out2.rmf")
#                                            self.get_tmp_file_name("out.rmf"));
        # verify that anchors remain intact during optimization
        if IMP.get_check_level() >= IMP.USAGE_AND_INTERNAL:
            short_init_factor=0.0001
            n_opt_cycles= 10
            print("short position initialization in non-fast mode")
        else:
            n_opt_cycles=100000
            short_init_factor=1.0
        IMP.npctransport.initialize_positions(sd,[],False,short_init_factor)

        print("Energy before optimization:",
              sd.get_bd().get_scoring_function().evaluate(False))

        print("Optimizing")
        print("dT %.1f fs" % sd.get_bd().get_maximum_time_step())
        sd.get_bd().optimize(n_opt_cycles)
        print()

        print("Energy after optimization:",
              sd.get_bd().get_scoring_function().evaluate(False))

        print("Report particles:")
        self._do_particles_report(sd.get_model(), sd.get_beads())
        print()

        sim_time_ns=500*TAU_NS
        sd.get_bd().set_current_time(0.0)
        sd.activate_statistics()
        print("Running for %.1f ns" % sim_time_ns)
        R=test_util.optimize_in_chunks(sd, sim_time_ns, NS_PER_CHUNK)
        #        print(R[0:100])

        self._analyze_run(assign_file, R)

    def test_harmonic_spring_score_from_protobuf_two_beads(self):
        '''
        test harmonic spring score and statistics for FG chains
        with two beads
        '''
        ntrials=3
        is_raise_on_fail= True
        if IMP.get_check_level() >= IMP.USAGE_AND_INTERNAL:
            print("Checks mode - limited run for speed")
            ntrials= 1
            is_raise_on_fail= False
        for i in range(ntrials):
            try:
                print("Try #", i)
                self._test_harmonic_spring_score_from_protobuf(2)
                return
            except AssertionError as e:
                if is_raise_on_fail:
                    print("EXCEPTION CAUGHT Try #", i)
                    print(e)
                    print("==\n\n")
                    if i+1==ntrials:
                        raise

    def test_harmonic_spring_score_from_protobuf_three_beads(self):
        '''
        test harmonic spring score and statistics for FG chains
        with three beads
        '''
        ntrials=3
        is_raise_on_fail= True
        if IMP.get_check_level() >= IMP.USAGE_AND_INTERNAL:
            print("Checks mode - limited run for speed")
            ntrials= 1
            is_raise_on_fail= False
        for i in range(ntrials):
            try:
                print("Try #", i)
                self._test_harmonic_spring_score_from_protobuf(3)
                return
            except AssertionError as e:
                if is_raise_on_fail:
                    print("EXCEPTION CAUGHT Try #", i)
                    print(e)
                    print("==\n\n")
                    if i+1==ntrials:
                        raise


if __name__ == '__main__':
    IMP.test.main()
