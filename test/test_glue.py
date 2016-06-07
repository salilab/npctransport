from __future__ import print_function
import IMP
import IMP.test
import IMP.npctransport
import RMF
import IMP.container
import math

# //r = rx * sqrt(s+1)
# r = rx^2 + ry^2
# s = ry^2/rx^2
# dU = -0.0625 * k * rx^2 * ry^2

radius=30
k_G=0.005
k=k_G
site_range=.5*radius
k_skew_G=0.25
rangeN = 10
rangeT = 20
range_skew_G = (rangeT/rangeN)**2
site_range_G=math.sqrt(rangeN**2 + rangeT**2)
kN = k_G*k_skew_G/(k_skew_G+1)
kT = k_G/(k_skew_G+1)
#fmax = max[(0.5 * kN * rangeN) * (0.25 * kT * rangeT ^ 2),
#           (0.5 * kT * rangeT) * (0.25 * kN * rangeN ^ 2)]
fmax = 0.125 * k_G * rangeN * rangeT * max(rangeN, rangeT)
dUmax = 0.0625 * k_G * rangeN**2 * rangeT**2
k=dUmax/site_range
nonspec_k_G=0.1*fmax
nonspec_range_G=.2*radius

print ("fmax=", fmax)
print ("dUmax=", dUmax)
print ("range=", site_range_G)
print ("rangeN=", rangeN)
print ("rangeT=", rangeT)
print ("kN=", kN)
print ("kT=", kT)
print ("k=", k_G)
print ("k_skew=",  k_skew_G)
print ("k_nonspec=",  nonspec_k_G)
print ("dU_nonspec=", nonspec_range_G*nonspec_k_G)

class Tests(IMP.test.TestCase):
    def _create_particle(self, m):
        p= IMP.Particle(m)
        d= IMP.core.XYZR.setup_particle(p)
        d.set_radius(radius)
        IMP.atom.Hierarchy.setup_particle(p)
        IMP.atom.Mass.setup_particle(p, 1.0)
        IMP.core.RigidBody.setup_particle(p, IMP.algebra.ReferenceFrame3D())
        d.set_coordinates_are_optimized(True)
        IMP.atom.RigidBodyDiffusion.setup_particle(p)
        return d
    def _test_one(self, site_range, site_k, nonspec_range, nonspec_k,
                  soft_sphere_k, dt, ntrial=0):
        print ("Test one:",
               "site-range", site_range, "[A]",
               "site-k", site_k, "[kCal/mol/A]",
               "nonspec-range", nonspec_range, "[A]",
               "nonspec-k", nonspec_k, "[kCal/mol/A]",
               "repulsive-k", soft_sphere_k, "[kCal/mol/A]",
               "dT", dt)
        print("dG-site:", site_range*site_k, "[kCal/mol]")
        nsteps = math.ceil(25E+6 / dt);
        if(IMP.get_check_level() >= IMP.USAGE):
            nsteps /= 250
        m= IMP.Model()
        m.set_log_level(IMP.PROGRESS)
        ds= [self._create_particle(m) for i in range(0,2)]
        dsi= [x.get_particle_index() for x in ds];
        ds[0].set_coordinates(IMP.algebra.Vector3D(0,0,0))
        ds[1].set_coordinates(IMP.algebra.Vector3D(2.5*radius,0,0))
        distance = IMP.algebra.get_distance(ds[0].get_coordinates(),
                                            ds[1].get_coordinates())
        print ("Initial distance", distance)
        types=[IMP.core.ParticleType(d.get_name()+" type") for d in ds]
        for d in zip(types, ds):
            IMP.core.Typed.setup_particle(d[1], d[0])
        sites=( [ IMP.algebra.Vector3D(radius, 0,0) ],
                [ IMP.algebra.Vector3D(-math.sqrt(1.0)*radius,math.sqrt(0.0)*radius,0)
#                  , IMP.algebra.Vector3D(-math.sqrt(0.8)*radius,math.sqrt(0.2)*radius,0)
                ] )
        ps= IMP.npctransport.SitesPairScore(site_range, site_k,
                                            nonspec_range, nonspec_k,
                                            soft_sphere_k,
                                            IMP.npctransport.get_spheres_from_vectors(sites[0], 0.0),
                                            IMP.npctransport.get_spheres_from_vectors(sites[1], 0.0))
#        ps.set_log_level(IMP.VERBOSE)
        r= IMP.core.PairRestraint(m, ps, dsi)
        bd= IMP.atom.BrownianDynamics(m)
        bd.set_scoring_function(r)
        bd.set_maximum_time_step(dt)
        f= RMF.create_rmf_file(self.get_tmp_file_name("glue1_%d.rmf" % ntrial))
        for d in zip(types, sites):
            IMP.npctransport.add_test_sites(f, d[0], site_range, d[1])
        IMP.rmf.add_restraint(f,r)
        w= IMP.npctransport.add_hierarchies_with_sites(f, ds)
        sos= IMP.rmf.SaveOptimizerState(m, f)
#        bd.add_optimizer_state(sos)
        sos.set_period(1000)
        max_delta = 0.75 * site_range
        sos.update_always()
        for rr,s in zip([2.00, 2.25, 2.50],
                       [-k*radius/2.0-0.1*k*0.2*radius,-k*radius/4.0, 0.0]):
            print ("rr*radius= ", rr*radius)
            ds[1].set_coordinates(IMP.algebra.Vector3D(rr*radius,0,0))
            print ("P0", ds[0])
            print ("P1", ds[1])
            init_score =  bd.get_scoring_function().evaluate(False)
            print("Initial score x0=", rr, "*R is ", init_score)
            self.assertAlmostEqual(init_score, s, delta = 0.001)
        for i in range(10):
            bd.optimize(math.ceil(nsteps / 10))
            sos.update_always()
            distance = IMP.algebra.get_distance(ds[0].get_coordinates(),
                                                ds[1].get_coordinates())
            print("score", i, " = ", bd.get_scoring_function().evaluate(False))
            if abs(distance - 2*radius) < max_delta:
                break;
        final_score = bd.get_scoring_function().evaluate(False)
        print("Final distance", distance, "score", final_score)
#        if(IMP.get_check_level() < IMP.USAGE):
        self.assertAlmostEqual(distance, 2*radius, delta =max_delta)
        self.assertLess(final_score, -0.001)
    def test_one(self):
        """Check interaction score repulsion for glue test"""
        IMP.set_log_level(IMP.PROGRESS)
        dt=IMP.npctransport.get_time_step(1, k, radius, .2*radius,0.1)
        ntrials=3
        for i in range(ntrials):
            try:
                self._test_one(site_range=.5*radius, site_k=k,
                            nonspec_range=.2*radius, nonspec_k=.1*k,
                               soft_sphere_k=0.1*k, dt=dt, ntrial=i)
                return
            except AssertionError:
                if i==(ntrials-1): raise

    def _test_two(self, site_range, site_k, nonspec_range, nonspec_k,
                  soft_sphere_k, dt):
        nsteps = 2000
        if(IMP.get_check_level() >= IMP.USAGE):
            nsteps /= 10
        m= IMP.Model()
        m.set_log_level(IMP.SILENT)
        ds= [self._create_particle(m) for i in range(0,3)]
        dsi= [x.get_particle_index() for x in ds];
        ds[0].set_coordinates(IMP.algebra.Vector3D(0,0,0))
        ds[1].set_coordinates(IMP.algebra.Vector3D(2*radius,0,0))
        ds[2].set_coordinates(IMP.algebra.Vector3D(0,0,2*radius))
        types=[IMP.core.ParticleType(d.get_name()+" type") for d in ds]
        for d in zip(types, ds):
            IMP.core.Typed.setup_particle(d[1], d[0])
        sites=([IMP.algebra.Vector3D(radius, 0,0),
                IMP.algebra.Vector3D(0, 0,radius)], [IMP.algebra.Vector3D(-radius, 0,0)],
               [IMP.algebra.Vector3D(0, 0,-radius)])
        rs=[]
        for p in [(0,1), (1,2), (0,2)]:
            ps= IMP.npctransport.SitesPairScore(site_range, site_k, nonspec_range,
                                                nonspec_k, soft_sphere_k,
                                                IMP.npctransport.get_spheres_from_vectors(sites[p[0]], 0.0),
                                                IMP.npctransport.get_spheres_from_vectors(sites[p[1]], 0.0))
            #          ps.set_log_level(IMP.VERBOSE)
            r=IMP.core.PairRestraint(m, ps, (dsi[p[0]], dsi[p[1]]))
            rs.append(r)
        bd= IMP.atom.BrownianDynamics(m)
        bd.set_scoring_function(rs)
        bd.set_maximum_time_step(dt)
        f= RMF.create_rmf_file(self.get_tmp_file_name("glue2.rmf"))
        for d in zip(types, sites):
            IMP.npctransport.add_test_sites(f, d[0], site_range, d[1])
        w= IMP.npctransport.add_hierarchies_with_sites(f, ds)
        IMP.rmf.add_restraints(f, rs)
        sos= IMP.rmf.SaveOptimizerState(m, f)
        sos.set_period(10)
        bd.add_optimizer_state(sos)
        bd.optimize(nsteps)
        sos.update_always()
    def test_two(self):
        """Check two interactions"""
        IMP.set_log_level(IMP.SILENT)
        dt=IMP.npctransport.get_time_step(1, k, radius,.2*radius)
        self._test_two(site_range=radius, site_k=k,
                       nonspec_range=.2*radius, nonspec_k=.5*k,
                       soft_sphere_k=k, dt=dt)


    def _create_restraint_three(self, m, ds, site_range, site_k, nonspec_range, nonspec_k,
                                soft_sphere_k, f=None):
        dsi= [x.get_particle_index() for x in ds];
        sites=([IMP.algebra.Vector3D(radius, 0,0), IMP.algebra.Vector3D(0, 0,radius)],
               [IMP.algebra.Vector3D(-radius, 0,0), IMP.algebra.Vector3D(0, 0,radius)],
               [IMP.algebra.Vector3D(radius, 0, 0), IMP.algebra.Vector3D(0, 0,-radius)],
               [IMP.algebra.Vector3D(-radius, 0,0), IMP.algebra.Vector3D(0, 0,-radius)])
        if f:
            for d in zip(ds, sites):
                IMP.npctransport.add_test_sites(f, IMP.core.Typed(d[0]).get_type(),
                                           .5*radius, d[1])
        rs=[]
        for p in [(0,1), (0,2), (0,3), (1,2), (1,3), (2,3)]:
            ps= IMP.npctransport.SitesPairScore(site_range, site_k, nonspec_range,
                                                nonspec_k, soft_sphere_k,
                                                IMP.npctransport.get_spheres_from_vectors(sites[p[0]], 0.0),
                                                IMP.npctransport.get_spheres_from_vectors(sites[p[1]], 0.0))
#          ps.set_log_level(IMP.VERBOSE)
            r=IMP.core.PairRestraint(m, ps, (dsi[p[0]], dsi[p[1]]))
            rs.append(r)
            r.evaluate(False)
            d=r.create_current_decomposition()
            if d:
                d.set_was_used(False)
        return rs
    def _test_three(self, site_range, site_k, nonspec_range, nonspec_k,
                  soft_sphere_k, dt):
        m= IMP.Model()
        m.set_log_level(IMP.SILENT)
        ds= [self._create_particle(m) for i in range(0,4)]
        dsi= [x.get_particle_index() for x in ds];
        ds[0].set_coordinates(IMP.algebra.Vector3D(0,0,0))
        ds[1].set_coordinates(IMP.algebra.Vector3D(2*radius,0,0))
        ds[2].set_coordinates(IMP.algebra.Vector3D(0,0,2*radius))
        ds[3].set_coordinates(IMP.algebra.Vector3D(2*radius,0,2*radius))
        types=[IMP.core.ParticleType(d.get_name()+" type") for d in ds]
        for d in zip(types, ds):
            IMP.core.Typed.setup_particle(d[1], d[0])
        f= RMF.create_rmf_file(self.get_tmp_file_name("glue3.rmf"))
        rs= self._create_restraint_three(m, ds,  site_range, site_k, nonspec_range, nonspec_k,
                                          soft_sphere_k, f)
        bd= IMP.atom.BrownianDynamics(m)
        bd.set_scoring_function(rs)
        bd.set_maximum_time_step(dt)
        w= IMP.npctransport.add_hierarchies_with_sites(f, ds)
        IMP.rmf.add_restraints(f, rs)
        sos= IMP.rmf.SaveOptimizerState(m, f)
        sos.set_period(1000)
        bd.add_optimizer_state(sos)
        print("optimizin")
        IMP.set_log_level(IMP.SILENT)
        bd.optimize(3000)
        print("done")
        sos.update_always()
    def test_three(self):
        """Check three interactions"""
        IMP.set_log_level(IMP.SILENT)
        dt=IMP.npctransport.get_time_step(1, k, radius, .2*radius)
        self._test_three(site_range=radius, site_k=k,
                       nonspec_range=.2*radius, nonspec_k=.5*k,
                       soft_sphere_k=k, dt=dt)
    def _rescore_three(self):
        IMP.set_log_level(IMP.SILENT)
        f= RMF.open_rmf_file_read_only(self.get_tmp_file_name("glue3.rmf"))
        m= IMP.Model()
        IMP.set_log_level(IMP.SILENT)
        ds= IMP.rmf.create_hierarchies(f, m)
        k=500
        rs= self._create_restraint_three(m, ds,  radius, k, .2*radius, .5*k, k)
        IMP.rmf.load_frame(f, RMF.FrameID(0))
        print(rs[0].evaluate(True))
        IMP.rmf.load_frame(f, RMF.FrameID(f.get_number_of_frames()-1))
        print (rs[0].evaluate(True))

    def _test_one_sliding(self, site_range, site_k,
                          range_skew, k_skew,
                          nonspec_range, nonspec_k,
                          soft_sphere_k, dt, ntrial=0):
        nsteps = 5E+6 / dt
        if(IMP.get_check_level() >= IMP.USAGE):
            nsteps /= 250
        m= IMP.Model()
        m.set_log_level(IMP.PROGRESS)
        ds= [self._create_particle(m) for i in range(0,2)]
        dsi= [x.get_particle_index() for x in ds];
        ds[0].set_coordinates(IMP.algebra.Vector3D(0,0,0))
        ds[1].set_coordinates(IMP.algebra.Vector3D(2.5*radius,0,0))
        ds[0].set_radius(radius*0.4) #/2.0*0.9) # for repulsion
        ds[1].set_radius(radius*0.8) # for repulsion
        distance = IMP.algebra.get_distance(ds[0].get_coordinates(),
                                            ds[1].get_coordinates())
        print ("Initial distance", distance)
        types=[IMP.core.ParticleType(d.get_name()+" type") for d in ds]
        for d in zip(types, ds):
            IMP.core.Typed.setup_particle(d[1], d[0])
        sites0 =  [ IMP.algebra.Sphere3D
                    ( IMP.algebra.Vector3D(0,0,0),
                      radius/2.0 ) ]
        sites1 = [ IMP.algebra.Sphere3D
                   ( IMP.algebra.Vector3D(-math.sqrt(1.0)*radius,
                                          math.sqrt(0.0)*radius,
                                          0.0),
                     0.0 ) ]
#                  , IMP.algebra.Vector3D(-math.sqrt(0.8)*radius,math.sqrt(0.2)*radius,0)
        print (site_range, site_k, range_skew, k_skew)
        print (nonspec_range, nonspec_k, soft_sphere_k)
        ps= IMP.npctransport.SitesPairScore(site_range, site_k,
                                            range_skew, k_skew,
                                            nonspec_range,
                                            nonspec_k,
                                            soft_sphere_k,
                                            sites0, sites1)
#        ps.set_log_level(IMP.VERBOSE)
        r= IMP.core.PairRestraint(m, ps, dsi)
        bd= IMP.atom.BrownianDynamics(m)
        bd.set_scoring_function(r)
        bd.set_maximum_time_step(dt)
        rmf_fname=self.get_tmp_file_name("glue1_slide_%d.rmf" % ntrial)
        rmf_fname="glue1_slide_%d.rmf" % ntrial
        print ("RMF:",rmf_fname)
        f= RMF.create_rmf_file(rmf_fname)
        for d in zip(types, [sites0, sites1]):
            IMP.npctransport.add_test_sites(f, d[0], d[1])
        IMP.rmf.add_restraint(f,r)
        w= IMP.npctransport.add_hierarchies_with_sites(f, ds)
        sos= IMP.rmf.SaveOptimizerState(m, f)
        bd.add_optimizer_state(sos)
        sos.set_period(math.ceil(100000.0/dt))
        max_delta = 0.75 * site_range
        sos.update_always()
        rangex = site_range_G / math.sqrt(range_skew_G+1)
        rangey = rangex * math.sqrt(range_skew_G)
        for rr,s in zip([1.50, 1.51, 1.8],
                       [-0.0625*k*(rangex*rangey)**2-nonspec_k*0.2*radius, 0.0, 0.0]):
            print ("rr*radius= ", rr*radius, " Expected score", s)
            ds[1].set_coordinates(IMP.algebra.Vector3D(rr*radius,0,0))
            print ("P0", ds[0])
            print ("P1", ds[1])
            init_score =  bd.get_scoring_function().evaluate(False)
            print ("Initial score [x0=", rr, "*R] is ", init_score)
#            self.assertAlmostEqual(init_score, s, delta = 0.001)
        ds[1].set_coordinates(IMP.algebra.Vector3D(1.4*radius,0,0))
        for i in range(10):
            bd.optimize(math.ceil(nsteps / 10))
            sos.update_always()
            distance = IMP.algebra.get_distance(ds[0].get_coordinates(),
                                                ds[1].get_coordinates())
            print ("score", i, " = ", bd.get_scoring_function().evaluate(False))
            if abs(distance - 2*radius) < max_delta:
                break;
        final_score = bd.get_scoring_function().evaluate(False)
        print ("Final distance", distance, "score", final_score)
#        if(IMP.get_check_level() < IMP.USAGE):
        self.assertAlmostEqual(distance, 2*radius, delta =max_delta)
        self.assertLess(final_score, -0.001)
    def test_one_sliding(self):
        """Check interaction score repulsion for glue test, sliding
        """
        print ("==\nSliding\n==")
        IMP.set_log_level(IMP.PROGRESS)
        dt=IMP.npctransport.get_time_step(1, k_G, radius, nonspec_range_G, 0.1)
        dt=min(dt,20000)
#        dt=50000
        print ("dT = ", dt)
        ntrials=3
        for i in range(ntrials):
            try:
                self._test_one_sliding(site_range=site_range_G,
                                       site_k=k_G,
                                       range_skew=range_skew_G,
                                       k_skew=k_skew_G,
                                       nonspec_range=nonspec_range_G,
                                       nonspec_k=nonspec_k_G,
                                       soft_sphere_k=100*math.sqrt(k_G/k_skew_G), dt=dt, ntrial=i)
                return
            except AssertionError:
                if i==(ntrials-1): raise



if __name__ == '__main__':
    IMP.test.main()
