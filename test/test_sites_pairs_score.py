from __future__ import print_function
import IMP
import IMP.test
import IMP.npctransport
import RMF
import IMP.container
import math
import numpy as np
from test_util import *

# r = rx * sqrt(s+1)
# r = rx^2 + ry^2
# s = ry^2/rx^2
# dU = -0.0625 * k * rx^2 * ry^2


radius=12
site_range=.5*radius
dUmax = 16.0 # kCal.mol
k_nonrot=dUmax/site_range; # set k_nonrot s.t. dUmax is the same for either the rotational or the non-rotational score
k_rot= 4*dUmax/site_range**2;
sigma1_max_deg=30
sigma2_max_deg=60
sigma1_max_rad=np.deg2rad(sigma1_max_deg)
sigma2_max_rad=np.deg2rad(sigma2_max_deg)
nonspec_k=0.01 #*k_nonrot
soft_sphere_k=10.0 #2*k_nonrot;
nonspec_range=site_range

print ("range=", site_range, "A")
print ("dUmax=", dUmax, " kCal/mol")
print ("k_rot=", k_rot, " kCal/mol/A^2")
print ("k_nonrot=", k_nonrot, " kCal/mol/A")
print ("soft_sphere_k=", soft_sphere_k, " kCal/mol/A")
print ("nonspec_k=", nonspec_k, " kCal/mol/A")
print ("nonspec_range", nonspec_range, "A")

def get_sliding_k_factor(s1,s2,s1max,s2max):
    ''' all in radians - estimate k factor in sites.h sliding score, independently of sites.h '''
    if(s1>s1max or s2>s2max):
        return 0
    k1=(math.cos(s1)-math.cos(s1max))/(1.0-math.cos(s1max))
    k2=(math.cos(s2)-math.cos(s2max))/(1.0-math.cos(s2max))
    return k1*k2

def print_info(ps, site0, site1):
    for p_serial,p in enumerate(ps):
        xyz=IMP.core.XYZ(p)
        xyzr=IMP.core.XYZR(p)
        rb=IMP.core.RigidBody(p)
        gForce=[xyz.get_derivative(ii) for ii in range(3)]
        lTorque=rb.get_torque()
        gTorque=rb.get_reference_frame().get_global_coordinates(lTorque)-xyz.get_coordinates()               
        if(p_serial==0):
            lSite=site0.get_center()
        else:
            lSite=site1.get_center()
        gSite=rb.get_reference_frame().get_global_coordinates(lSite)
        print (p.get_name(), "gCoords", xyz,"gSite",gSite,"gForce",gForce,"gTorque",gTorque)

class Tests(IMP.test.TestCase):
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
        nsteps = math.ceil(50E+6 / dt);
        if(IMP.get_check_level() >= IMP.USAGE):
            nsteps /= 250
        m= IMP.Model()
        ps= [create_diffusing_rb_particle(m,radius) for i in range(0,2)]
        ds= [IMP.core.XYZR(p) for p in ps];
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
                                            0.0, 0.0,
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
        for x_ds1,estimated_score in [
            (2*radius, -k_nonrot*site_range - nonspec_k*nonspec_range),
            (2*radius+0.5*site_range,  -0.5*k_nonrot*site_range - 0.5*nonspec_k*nonspec_range),
            (-2*radius,-nonspec_k*nonspec_range),
            (radius, soft_sphere_k*radius-nonspec_k*nonspec_range),
            (2*radius+site_range, 0.0)
            ]:
            ds[1].set_coordinates(IMP.algebra.Vector3D(x_ds1,0,0))
#            print ("P0", ds[0])
#            print ("P1", ds[1])
            init_score =  bd.get_scoring_function().evaluate(False)
            print("Score for sphere-distance ", abs(x_ds1)-2*radius, 
                  "A and site-distance", abs(x_ds1-2*radius),
                  "A is: %.2f vs. expected: %.2f" % (init_score, estimated_score)
                  )
            self.assertAlmostEqual((init_score+0.0001)/(estimated_score+0.0001), 1.0, delta = 0.0001)
        print("Score before optimization is %.2f" % init_score)
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
        IMP.set_log_level(IMP.SILENT)
        dt=IMP.npctransport.get_time_step(1, 
                                          k_nonrot, 
                                          radius, 
                                          nonspec_range,  
                                          0.025)
        print("TEST ONE dt=", dt)
        ntrials=3
        for i in range(ntrials):
            try:
                self._test_one(site_range=site_range, site_k=k_nonrot,
                            nonspec_range=nonspec_range, nonspec_k=nonspec_k,
                               soft_sphere_k=soft_sphere_k, dt=dt, ntrial=i)
                return
            except AssertionError:
                if i==(ntrials-1): raise

    def _test_two(self, site_range, site_k, nonspec_range, nonspec_k,
                  soft_sphere_k, dt):
        nsteps = 2000
        if(IMP.get_check_level() >= IMP.USAGE):
            nsteps /= 10
        m= IMP.Model()
        ps= [create_diffusing_rb_particle(m,radius) for i in range(0,3)]
        ds= [IMP.core.XYZR(p) for p in ps];
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
            ps= IMP.npctransport.SitesPairScore(site_range, site_k, 
                                                0.0, 0.0,
                                                nonspec_range, nonspec_k, soft_sphere_k,
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
        dt=IMP.npctransport.get_time_step(1, 
                                          k_nonrot, 
                                          radius, 
                                          nonspec_range,  
                                          0.1)
        self._test_two(site_range=site_range, site_k=k_nonrot,
                       nonspec_range=nonspec_range, nonspec_k=nonspec_k,
                       soft_sphere_k=soft_sphere_k, dt=dt)


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
            ps= IMP.npctransport.SitesPairScore(site_range, site_k, 
                                                0.0, 0.0,
                                                nonspec_range, nonspec_k, soft_sphere_k,
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
        ps= [create_diffusing_rb_particle(m,radius) for i in range(0,4)]
        ds= [IMP.core.XYZR(p) for p in ps];
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
        dt=IMP.npctransport.get_time_step(1,
                                          k_nonrot,
                                          radius,
                                          nonspec_range,
                                          0.1)
        self._test_three(site_range=site_range, site_k=k_nonrot,
                         nonspec_range=nonspec_range, nonspec_k=nonspec_k,
                         soft_sphere_k=soft_sphere_k, dt=dt)
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
                          sigma1_max_deg, sigma2_max_deg,
                          nonspec_range, nonspec_k,
                          soft_sphere_k, dt, ntrial=0):
        time_ns=50
        nsteps = time_ns * 1E+6 / dt
        if(IMP.get_check_level() >= IMP.USAGE):
            nsteps /= 250
            print("FAST SLIDING")
        m= IMP.Model()
        ps= [create_diffusing_rb_particle(m,r)  \
             for r in [radius,0.5*radius] ]
        pis= [p.get_index() for p in ps];
        xyzrs= [IMP.core.XYZR(p) for p in ps];
        rbs= [IMP.core.RigidBody(p) for p in ps];
        xyzrs[0].set_coordinates(IMP.algebra.Vector3D(0,0,0))
        xyzrs[1].set_coordinates(IMP.algebra.Vector3D(2.5*radius,0,0))
        distance = IMP.algebra.get_distance(xyzrs[0].get_coordinates(),
                                            xyzrs[1].get_coordinates())
        print ("Initial distance", distance)
        types=[IMP.core.ParticleType(p.get_name()+" type") for p in ps]
        for p, type in zip(ps, types):
            IMP.core.Typed.setup_particle(p,type)
        sites0 =  [ IMP.algebra.Sphere3D
                    ( IMP.algebra.Vector3D(radius,0,0),
                      0.0 ) ]
        sites1 = [ IMP.algebra.Sphere3D
                   ( IMP.algebra.Vector3D(-math.sqrt(1.0)*0.5*radius,
                                          math.sqrt(0.0)*0.5*radius,
                                          0.0),
                     0.0 ) ]
        print('Sites:', sites0, sites1)
#                  , IMP.algebra.Vector3D(-math.sqrt(0.8)*radius,math.sqrt(0.2)*radius,0)
        sps= IMP.npctransport.SitesPairScore(site_range, site_k,
                                            sigma1_max_deg, sigma2_max_deg,
                                            nonspec_range,
                                            nonspec_k,
                                            soft_sphere_k,
                                            sites0, sites1)
#        ps.set_log_level(IMP.VERBOSE)
        r= IMP.core.PairRestraint(m, sps, pis)
        bd= IMP.atom.BrownianDynamics(m)
        bd.set_scoring_function(r)
        bd.set_maximum_time_step(dt)
        rmf_fname=self.get_tmp_file_name("glue1_slide_%d.rmf" % ntrial)
        rmf_fname="glue1_slide_%d.rmf" % ntrial
#        print ("RMF:",rmf_fname)
        f= RMF.create_rmf_file(rmf_fname)
        for type,site in zip(types, [sites0, sites1]):
            print("Type/site", type,site)
            IMP.npctransport.add_test_sites(f, type, site)
        IMP.rmf.add_restraint(f,r)
        w= IMP.npctransport.add_hierarchies_with_sites(f, xyzrs)
        sos= IMP.rmf.SaveOptimizerState(m, f)
        bd.add_optimizer_state(sos)
        max_delta = 0.75 * site_range
        sos.update_always()
        print("TMP", np.rad2deg(math.acos(0.5+0.5*math.cos(sigma1_max_rad))))
        for i,p in enumerate(ps):
            lForce=[IMP.core.XYZ(p).get_derivative(ii) for ii in range(3)]
            lTorque=rbs[i].get_torque()
            print ("P%d" % i, "gCoords", xyzrs[i],"lForce",lForce,"lTorque",lTorque)
        for rr,s,sigma0,sigma1 in [
            (1.50,
             -0.25*k_rot*site_range**2 - nonspec_k*nonspec_range,
             0.0,
             0.0),
            (1.50,
             -0.5*0.25*k_rot*site_range**2 - nonspec_k*nonspec_range,
             math.acos(0.5+0.5*math.cos(sigma1_max_rad)),
             0.0),
            (1.50,
             -0.5*0.25*k_rot*site_range**2 - nonspec_k*nonspec_range,
             0.0,
             math.acos(0.5+0.5*math.cos(sigma2_max_rad))),
            (1.50, -0.25*0.25*k_rot*site_range**2-nonspec_k*nonspec_range,
             math.acos(0.5+0.5*math.cos(sigma1_max_rad)),
             math.acos(0.5+0.5*math.cos(sigma2_max_rad))),
             (1.50+0.5*site_range/radius,
              -0.125*k_rot*site_range**2 - 0.5*nonspec_k*nonspec_range,
              0.0,
              0.0),
             (1.50+0.5*site_range/radius,
              -get_sliding_k_factor(sigma1_max_rad/2.0, sigma2_max_rad/3.0,sigma1_max_rad,sigma2_max_rad)*0.125*k_rot*site_range**2 - 0.5*nonspec_k*nonspec_range,
              sigma1_max_rad/2.0,
              sigma2_max_rad/3.0)
             ]:
            print ("--\nrr*radius= ", rr*radius,
                   "sigmas", sigma0/3.141*180, sigma1/3.141*180,
                   "Expected score", s)
            rot_vector=IMP.algebra.Vector3D(0,0,1)
            tr1=rbs[1].get_reference_frame().get_transformation_to()
            rot0=IMP.algebra.get_rotation_about_normalized_axis(rot_vector,sigma0)
            rot1=IMP.algebra.get_rotation_about_normalized_axis(rot_vector,sigma1)
            tr0=IMP.algebra.Transformation3D(rot0,IMP.algebra.Vector3D(0,0,0))
            tr1=IMP.algebra.Transformation3D(rot1,IMP.algebra.Vector3D((rr*radius,0,0)))
            rbs[0].set_reference_frame(IMP.algebra.ReferenceFrame3D(tr0));
            rbs[1].set_reference_frame(IMP.algebra.ReferenceFrame3D(tr1));
            print("Reference Frames:")
            print(rbs[0].get_reference_frame())
            print(rbs[1].get_reference_frame())
            print("1,0,0 transformed by rb0 ref-frame: ",
                  rbs[0].get_reference_frame().get_global_coordinates(IMP.algebra.Vector3D(1,0,0)))
            print(rbs[0].get_reference_frame().get_transformation_to().get_rotation().get_rotated(IMP.algebra.Vector3D(1,0,0)))
            init_score =  bd.get_scoring_function().evaluate(True)
            print ("Computed score [x0=", rr, "*R] is ", init_score)
            if(s is not None):
                self.assertAlmostEqual(init_score, s, delta = 0.001)
            print_info(ps,sites0[0],sites1[0])
            sos.update_always()
    
        xyzrs[1].set_coordinates(IMP.algebra.Vector3D(1.5*radius+0.5*site_range,0,0))
        sos.set_period(math.ceil(nsteps/100000))
        for i in range(10):
            print("Steps",nsteps)
            bd.optimize(math.ceil(nsteps / 10))
            sos.update_always()
            distance = IMP.algebra.get_distance(xyzrs[0].get_coordinates(),
                                                xyzrs[1].get_coordinates())
            print ("score", i, " = ", bd.get_scoring_function().evaluate(False))
            if abs(distance - 1.5*radius) < max_delta:
                break;
        final_score = bd.get_scoring_function().evaluate(False)
        print ("Final distance", distance, "score", final_score)

        self.assertAlmostEqual(distance, 1.5*radius, delta =max_delta)
        self.assertLess(final_score, -0.001)

        
    def test_one_sliding(self):
        """Check interaction score repulsion for glue test, sliding
        """
        print ("==\nSliding\n==")
        IMP.set_log_level(IMP.SILENT)
        dt=IMP.npctransport.get_time_step(1,
                                          k_nonrot,
                                          radius,
                                          min(site_range,nonspec_range),
                                          0.05)
        dt=min(dt,2000)
#        dt=50000
        print ("dT = ", dt)
        ntrials=3
        for i in range(ntrials):
            try:
                self._test_one_sliding(site_range=site_range,
                                       site_k=k_rot,
                                       sigma1_max_deg=sigma1_max_deg,
                                       sigma2_max_deg=sigma2_max_deg,
                                       nonspec_range=nonspec_range,
                                       nonspec_k=nonspec_k,
                                       soft_sphere_k=soft_sphere_k, 
                                       dt=dt, ntrial=i)
                return
            except AssertionError:
                if i==(ntrials-1): raise



if __name__ == '__main__':
    IMP.test.main()
