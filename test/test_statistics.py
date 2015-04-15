from __future__ import print_function
import IMP
import IMP.test
import IMP.npctransport
import IMP.container
import math

radius=5

class Tests(IMP.test.TestCase):
    def _create_diffusing_particle(self, m, radius= 10.0):
        """creates a particle that is decorated with RigidBody,
        RigidBodyDiffusion and XYZR coordinates

        m - model
        """
        p= IMP.Particle(m)
        d= IMP.core.XYZR.setup_particle(p)
        d.set_radius( radius )
        IMP.core.RigidBody.setup_particle(p, IMP.algebra.ReferenceFrame3D())
        dd= IMP.atom.RigidBodyDiffusion.setup_particle(p)
        return p

    def test_repulsion(self):
        """Check diffusion coefficient estimation"""
        print("TEST_REPULSION")
        if IMP.base.get_check_level() >= IMP.base.USAGE_AND_INTERNAL:
            print("INTERNAL")
            n_cycles = 1000
            delta_factor = 0.5
        else:
            print("NORMAL")
            n_cycles=500000
            delta_factor = 0.1
        m= IMP.Model()
        p= self._create_diffusing_particle(m)
        IMP.core.XYZR(p).set_coordinates_are_optimized(True)
        dt=1000
        bd= IMP.atom.BrownianDynamics(m)
        bd.set_maximum_time_step(dt)
        os= IMP.npctransport.BodyStatisticsOptimizerState(p)
        os.set_period(10)
        bd.add_optimizer_state(os)
        IMP.set_log_level(IMP.SILENT)
        bd.set_log_level(IMP.SILENT)
        bd.optimize(n_cycles)
        Dout= os.get_diffusion_coefficient()
        Din= IMP.atom.RigidBodyDiffusion(p).get_diffusion_coefficient()
        print(Dout, Din)
        self.assertAlmostEqual(Dout, Din,
                               delta=delta_factor*Din)
    def test_rot(self):
        """Check rigid body correlation time"""
        print("TEST_ROT")
        if IMP.build!= "fast":
          self.skipTest("Only run in fast mode")
        m= IMP.Model()
        p= self._create_diffusing_particle(m)
        IMP.core.RigidBody(p).set_coordinates_are_optimized(True)
        dd= IMP.atom.RigidBodyDiffusion(p)
        nD=dd.get_rotational_diffusion_coefficient()
        dd.set_rotational_diffusion_coefficient(10*nD)
        print(dd.get_rotational_diffusion_coefficient() \
            , dd.get_diffusion_coefficient())
        dt=100000
        bd= IMP.atom.BrownianDynamics(m)
        bd.set_maximum_time_step(dt)
        os= IMP.npctransport.BodyStatisticsOptimizerState(p)
        num_steps=1000
        os.set_period(num_steps/1000)
        bd.add_optimizer_state(os)
        IMP.set_log_level(IMP.SILENT)
        bd.set_log_level(IMP.SILENT)
        bd.optimize(num_steps)
        cor_out= os.get_correlation_time()
        Din= dd.get_rotational_diffusion_coefficient()
        Dout=1.0/(2.0*cor_out)
        print(Dout, Din, cor_out)
        self.assertAlmostEqual(Dout, Din, delta=.5*Dout)
    def _create_magnet_restraint(self, m, p, magnet_coordinates):
        """
        create a dummy particle at magnet_coordinate,
        and return a restraint that magnetizes p towards it

        m - the model
        p - the particle to be magnetized towards magnet_coordinates
        magnet_coordinates - the coordinates of the dummy magnet particle
        """
        p_magnet= IMP.Particle(m)
        p_magnet_rb= IMP.core.XYZ.setup_particle(p_magnet, magnet_coordinates)
        p_magnet_rb.set_coordinates_are_optimized(False)
        harmonic = IMP.core.Harmonic(0,1)
        magnet_restraint= IMP.core.DistanceRestraint(harmonic, p, p_magnet)
        return (magnet_restraint, p_magnet_rb)
    def test_transport_stats(self):
        """Check particle transport stats"""
        print("TEST_TRANSPORT_STATS")
        print("\nTesting particle transport statistics:")
        if IMP.build!= "fast":
          self.skipTest("Only run in fast mode")
        m= IMP.Model()
        p= self._create_diffusing_particle(m)
        p_rb= IMP.core.RigidBody(p)
        p_rb.set_coordinates([0,0,0])
        p_rb.set_coordinates_are_optimized(True)
        # create a magnet that pulls p towards it:
        (magnet_restraint, p_magnet_rb)= self._create_magnet_restraint(m, p, [0,0,30])
        m.add_restraint( magnet_restraint)
        # make simulation of p running towards p_magnet:
        dt= 10
        bd= IMP.atom.BrownianDynamics(m)
        bd.set_maximum_time_step(dt)
        print("hey")
        os= IMP.npctransport.ParticleTransportStatisticsOptimizerState(p,10,20)
        print("ho")
        num_steps=10000
        os.set_period( num_steps / 1000 )
        bd.add_optimizer_state( os )
        IMP.set_log_level(IMP.SILENT)

        print("Before optimization z = %.2f" % p_rb.get_coordinates()[2])
        bd.optimize( num_steps )
        print("After 1st optimization z = %.2f" % p_rb.get_coordinates()[2])
        p_magnet_rb.set_coordinates( [0,0,0] );
        bd.optimize( num_steps )
        print("After 2nd optimization z = %.2f" % p_rb.get_coordinates()[2])
        n = os.get_total_n_transports()
        print("# of transports: %d" % n)
        self.assertEqual(n,2)
    def test_rot_nrb(self):
        """Check hidden rigid body correlation time"""
        print("TEST_ROT_NRB")
        m= IMP.Model()
        p= IMP.Particle(m, "rb")
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
        cr= IMP._ConstRestraint(0, rb.get_rigid_members())
        m.add_restraint(cr)
        dd= IMP.atom.RigidBodyDiffusion.setup_particle(p)
        nD=10.0*dd.get_rotational_diffusion_coefficient()
        dd.set_rotational_diffusion_coefficient(nD)
        dt=10000
        bd= IMP.atom.BrownianDynamics(m)
        bd.set_maximum_time_step(dt)
        os= IMP.npctransport.ChainStatisticsOptimizerState(rb.get_rigid_members())
        os.set_period(10)
        bd.add_optimizer_state(os)
        os2= IMP.npctransport.BodyStatisticsOptimizerState(p)
        os2.set_period(10)
        bd.add_optimizer_state(os2)
        IMP.set_log_level(IMP.SILENT)
        bd.set_log_level(IMP.SILENT)
        bd.optimize(5000)
        Dout= os.get_correlation_time()
        Dout2= os2.get_correlation_time()
        Din= dd.get_rotational_diffusion_coefficient()
        v=1.0/(2.0*Dout)
        print("Corr1, Corr2, rot_dif_coeff, v",  Dout, Dout2, Din, v)
        self.assertAlmostEqual(Dout, Dout2, delta=.1*Dout)
        dfs= os.get_diffusion_coefficients()
        print(dfs)
        for d in dfs:
            self.assertAlmostEqual(0, d, delta=.1)

if __name__ == '__main__':
    IMP.test.main()
