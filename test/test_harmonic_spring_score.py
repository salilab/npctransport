from __future__ import print_function
import IMP
import IMP.test
import IMP.npctransport
import IMP.container
import IMP.rmf
import RMF
import math
import numpy
import time
import random
import test_util

radius=7

def do_particles_report(m, pis):
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


class ConeTests(IMP.test.TestCase):
    def test_harmonic_spring_score(self):
        ntrials=3
        for i in range(ntrials):
            try:
                print("Try #", i)
                self._test_harmonic_spring_score()
                return
            except:
                print("EXCEPTION CAUGHT Try #", i)
                f= RMF.create_rmf_file(self.get_tmp_file_name("tmp%d.rmf" % i)) # just to force a flush
                if i+1==ntrials:
                    raise
    def _create_diffuser(self, m):
        p= IMP.Particle(m)
        d= IMP.core.XYZR.setup_particle(p)
        d.set_radius(radius)
        d.set_coordinates_are_optimized(True)
        h= IMP.atom.Hierarchy.setup_particle(p)
        m= IMP.atom.Mass.setup_particle(p, 1)
        diff= IMP.atom.Diffusion.setup_particle(p)
        return d
    def _randomize(self, ds, bb):
        for d in ds:
            d.set_coordinates(IMP.algebra.get_random_vector_in(bb))
    def _show(self, ds, w):
        for d in ds:
            g= IMP.core.XYZRGeometry(d);
            w.add_geometry(g)
    def _test_harmonic_spring_score(self):
        """Check linear well"""
        m= IMP.Model()
        T=298
        m.set_log_level(IMP.SILENT)
        ds= [self._create_diffuser(m) for i in range(0,2)]
        pis=[x.get_particle_index() for x in ds]
        ds[1].set_coordinates(IMP.algebra.Vector3D(0,2*radius,0))
        rest_length_factor = 1.75
        k = 0.5 # kcal/mol/A^2
        rest_length=rest_length_factor*radius*2.0
        tau_ns= 1
        tau_fs= tau_ns*(1E+6)
        D_A2_per_fs= IMP.atom.get_kt(T)/tau_fs/k
        print ("D [A^2/fs]", D_A2_per_fs, "D particle", IMP.atom.Diffusion(m,pis[0]).get_diffusion_coefficient())
        rs= IMP.npctransport.RelaxingSpring.setup_particle \
            (m,
             pis[0],
             pis[0],
             pis[1],
             rest_length,
             D_A2_per_fs)
        ss= IMP.npctransport.HarmonicSpringSingletonScore(20*k, k)
        r= IMP.core.SingletonRestraint(m, ss, pis[0])
        bd= IMP.npctransport.BrownianDynamicsTAMDWithSlabSupport(m)
        bd.set_maximum_time_step(1000)
        bd.set_scoring_function([r])
        bd.set_temperature(T)
        rmf_fname= self.get_tmp_file_name("well%d.rmf" % int(random.random()*1000))
        rmf_fname= "test.rmf"
        print(rmf_fname)
        f= RMF.create_rmf_file(rmf_fname)
        # TODO: note that if ds would have contained sites, we would
        # need to switch to npctransport::add_hierarchies_with_sites(),
        # perhaps worth switching anyway?
        IMP.rmf.add_hierarchies(f, ds)
        IMP.rmf.add_restraints(f, [r])
        w= IMP.rmf.SaveOptimizerState(m, f)
        w.set_period(100) # set this to one only for debugging
        bd.add_optimizer_state(w)
        Edist=0.0
        Erest=0.0
        n=0.0
        inner=50*tau_ns
        if IMP.get_check_level() >= IMP.USAGE_AND_INTERNAL:
            outer=10
        else:
            outer=10000
        T=[0]
        D=[IMP.core.get_distance(ds[0], ds[1])]
        R=[rs.get_rest_length()]
        do_particles_report(m, pis)
        for i in range(outer):
#            print("%.1f [ns]\t" % (bd.get_current_time()*1E-6),end='')
            bd.optimize(inner)
            dist= IMP.core.get_distance(ds[0], ds[1])
            rest_length= rs.get_rest_length()
            T.append(bd.get_current_time()*1E-6)
            D.append(dist)
            R.append(rest_length)
            n0= n
            n= n+1
            Edist= (Edist*n0+dist)/n
            Erest= (Erest*n0+rest_length)/n
#            print(dist, " ", rest_length)
        print("Edist %.2f" % Edist, "Erest %.2f" % Erest, "Eq-rest %.2f" % rs.get_equilibrium_rest_length())
        ExpectedDelta= 0.5 + math.sqrt(3/k/(bd.get_current_time()/tau_fs)) # delta scales with sqrt(3/k) for k in kcal/mol/A^2 beause the spring energy 0.5*k*R^2 is in the order of [kB]T, or ~0.6 kcal/mol, so 0.3*k*[kB]T should be on the order of [kB]T. The mean will converge with simulation time, but the order of [kB]T is a safe margin
        print("ExpectedDelta", ExpectedDelta)
        self.assertAlmostEqual(Erest,
                               rs.get_equilibrium_rest_length(),
                               delta=ExpectedDelta)

        if IMP.get_check_level() >= IMP.USAGE_AND_INTERNAL:
            return

        # Check that relaxation time is indeed on the order of tau
        # = autocorrelation decays exponentially with time/tau
        try:
            import pandas
        except ImportError:
            print("WARNING: pandas module not installed, skipping autocorrelation test")
            return
        Rp= pandas.Series(R)
        dT_fs= inner*bd.get_maximum_time_step()
        i= int(round(tau_fs/dT_fs))
        print("i", i, "tau_fs/dT_fs", tau_fs/dT_fs, "tau_ns", tau_ns, "dT_ns", dT_fs*1E-6)
        print("corr at tau_ns: %.3f" % Rp.autocorr(i))
        self.assertAlmostEqual(Rp.autocorr(i),
                               math.exp(-1),
                               delta=0.04)
#        for i in range(len(R)):
#            print (T[i],Rp.autocorr(i))



if __name__ == '__main__':
    IMP.test.main()
