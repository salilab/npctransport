import IMP
import IMP.algebra
import IMP.atom
import IMP.container
import IMP.core
import IMP.npctransport
import IMP.rmf
import IMP.test
import RMF
import random


class Tests(IMP.test.TestCase):

    """Test the GranuleActivationOptimizerState class"""

    def _create_particle(self, m, radius):
        p= IMP.Particle(m, "body")
        xyzr= IMP.core.XYZR.setup_particle(p)
        xyzr.set_radius(radius)
        xyzr.set_coordinates_are_optimized(True)
        rb= IMP.core.RigidBody.setup_particle(p, IMP.algebra.ReferenceFrame3D())
        rbd= IMP.atom.RigidBodyDiffusion.setup_particle(p)
#        IMP.atom.Mass.setup_particle(p, 1)
#        IMP.atom.Diffusion.setup_particle(p)
        IMP.display.Colored.setup_particle(p, IMP.display.get_display_color(2))
        return p


    def test_on_particles(self):
        """Check with Brownian dynamics of particles"""
        m= IMP.Model()
#        m.set_log_level(IMP.SILENT)
        ps0= []
        ps1= []
        n0= 10
        n1= 50
        for i in range(n0):
            ps0.append( self._create_particle(m, radius=10.0) )
        for i in range(n1):
            ps1.append( self._create_particle(m, radius=5.0) )

        rs=[] # restraints

        ps_excluded= IMP.core.SoftSpherePairScore(10.0)
        cpc= IMP.container.ClosePairContainer(ps0+ps1, 0.0, 3.0)
        rs.append( IMP.container.PairsRestraint(ps_excluded, cpc) )

        L= 400
        bb= IMP.algebra.BoundingBox3D(IMP.algebra.Vector3D(-L/2, -L/2, -L/2), IMP.algebra.Vector3D(L/2, L/2, L/2))
        bbss= IMP.core.BoundingBox3DSingletonScore(IMP.core.HarmonicUpperBound(0, 10),
                                                          bb)
        rs.append(IMP.container.SingletonsRestraint(bbss, ps0+ps1))

        sf= IMP.core.RestraintsScoringFunction(rs)

        for p in (ps0+ps1):
#            rot= IMP.algebra.get_random_rotation_3d()
            tr= IMP.algebra.get_random_vector_in(bb)
#            rf= IMP.algebra.ReferenceFrame3D(IMP.algebra.Transformation3D(rot,
#                                                                           tr))
#            IMP.core.RigidBody(pa).set_reference_frame(rf)
            IMP.core.XYZ(p).set_coordinates(tr)


        f= RMF.create_rmf_file("test.rmf")
        p_root= IMP.Particle(m, "root")
        h_root= IMP.atom.Hierarchy.setup_particle(p_root)
        for p in (ps0+ps1):
            h_child= IMP.atom.Hierarchy.setup_particle(p)
            h_root.add_child(h_child)
        IMP.rmf.add_hierarchy(f, h_root)
        IMP.rmf.add_restraints(f, rs)
        sos= IMP.rmf.SaveOptimizerState(m, f)
        sos.set_period(100)

        gaos= IMP.npctransport.GranuleActivationOptimizerState(ps0, ps1,
                                                               25,
                                                               5,
                                                               10)
        print("Model", gaos.get_model())
        IMP.set_log_level(IMP.SILENT)
#        gaos.set_log_level(IMP.TERSE)

        bd= IMP.atom.BrownianDynamics(m)
#        sf.set_log_level(IMP.SILENT)
        bd.set_scoring_function(sf)
        bd.set_maximum_time_step(100)
        bd.add_optimizer_state(sos)
        bd.add_optimizer_state(gaos)

#        bd.optimize(10)
#        print("going silent")
#        IMP.set_log_level(IMP.V)
        max_cycles= 5000000
        round_cycles= 50000
        total_cycles= 0
        e_threshold= 0.3
        for i in range(max_cycles / round_cycles):
            bd.optimize(round_cycles)
#            gaos.update()
            energy= sf.evaluate(False)
            total_cycles += round_cycles
            print("energy after %d cycles = %.2f" \
                % (total_cycles, energy))
            print("Bounding box score:")
            print(bbss.evaluate_indexes(m, ps0+ps1, None, 0, n0+n1-1))
            print("Excluded score:")
            print(ps_excluded.evaluate_indexes(m, cpc.get_indexes(), None,
                                                  0, len(cpc.get_indexes())))
            print("--")

if __name__ == '__main__':
    IMP.test.main()
