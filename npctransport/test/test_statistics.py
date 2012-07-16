import IMP
import IMP.test
import IMP.npctransport
import IMP.container
import math

radius=5

class Tests(IMP.test.TestCase):
    def test_repulsion(self):
        """Check diffusion coefficient estimation"""
        m= IMP.Model()
        p =IMP.Particle(m)
        d=IMP.core.XYZR.setup_particle(p)
        d.set_radius(10)
        d.set_coordinates_are_optimized(True)
        IMP.core.RigidBody.setup_particle(p, IMP.algebra.ReferenceFrame3D())
        dd= IMP.atom.RigidBodyDiffusion.setup_particle(p)
        dt=1000
        bd= IMP.atom.BrownianDynamics(m)
        bd.set_maximum_time_step(dt)
        os= IMP.npctransport.BodyStatisticsOptimizerState(p)
        os.set_period(10)
        bd.add_optimizer_state(os)
        bd.optimize(1000)
        Dout= os.get_diffusion_coefficient()
        Din= dd.get_diffusion_coefficient()
        print Dout, Din
        self.assertAlmostEqual(Dout, Din, delta=.4*Dout)
    def test_rot(self):
        """Check rigid body correlation time"""
        if IMP.build!= "fast":
          self.skipTest("Only run in fast mode")
        m= IMP.Model()
        p =IMP.Particle(m)
        d=IMP.core.XYZR.setup_particle(p)
        d.set_radius(10)
        rb= IMP.core.RigidBody.setup_particle(p, IMP.algebra.ReferenceFrame3D())
        dd= IMP.atom.RigidBodyDiffusion.setup_particle(p)
        rb.set_coordinates_are_optimized(True)
        nD=dd.get_rotational_diffusion_coefficient()
        dd.set_rotational_diffusion_coefficient(10*nD)
        print dd.get_rotational_diffusion_coefficient(), dd.get_diffusion_coefficient()
        dt=100000
        bd= IMP.atom.BrownianDynamics(m)
        bd.set_maximum_time_step(dt)
        os= IMP.npctransport.BodyStatisticsOptimizerState(p)
        num_steps=1000
        os.set_period(num_steps/1000)
        bd.add_optimizer_state(os)
        IMP.set_log_level(IMP.SILENT)
        bd.optimize(num_steps)
        cor_out= os.get_correlation_time()
        Din= dd.get_rotational_diffusion_coefficient()
        Dout=1.0/(2.0*cor_out)
        print Dout, Din, cor_out
        self.assertAlmostEqual(Dout, Din, delta=.5*Dout)

    def test_rot_nrb(self):
        """Check hidden rigid body correlation time"""
        m= IMP.Model()
        p =IMP.Particle(m, "rb")
        d=IMP.core.XYZR.setup_particle(p)
        d.set_radius(10)
        d.set_coordinates_are_optimized(True)
        rb= IMP.core.RigidBody.setup_particle(p, IMP.algebra.ReferenceFrame3D())
        for i in range(3):
            pc= IMP.Particle(m, str(i))
            d=IMP.core.XYZR.setup_particle(pc)
            d.set_radius(2)
            v= IMP.algebra.Vector3D(0,0,0)
            v[i]=5
            d.set_coordinates(v)
            rb.add_member(pc)
        # to make sure coordinates get updated
        cr= IMP._ConstRestraint(0, rb.get_members())
        m.add_restraint(cr)
        dd= IMP.atom.RigidBodyDiffusion.setup_particle(p)
        nD=10.0*dd.get_rotational_diffusion_coefficient()
        dd.set_rotational_diffusion_coefficient(nD)
        dt=10000
        bd= IMP.atom.BrownianDynamics(m)
        bd.set_maximum_time_step(dt)
        os= IMP.npctransport.ChainStatisticsOptimizerState(rb.get_members())
        os.set_period(10)
        bd.add_optimizer_state(os)
        os2= IMP.npctransport.BodyStatisticsOptimizerState(p)
        os2.set_period(10)
        bd.add_optimizer_state(os2)
        IMP.set_log_level(IMP.SILENT)
        bd.optimize(1000)
        Dout= os.get_correlation_time()
        Dout2= os2.get_correlation_time()
        Din= dd.get_rotational_diffusion_coefficient()
        v=1.0/(2.0*Dout)
        print Dout, Dout2, Din, v
        self.assertAlmostEqual(Dout, Dout2, delta=.1*Dout)
        dfs= os.get_diffusion_coefficients()
        print dfs
        for d in dfs:
            self.assertAlmostEqual(0, d, delta=.1)

if __name__ == '__main__':
    IMP.test.main()
