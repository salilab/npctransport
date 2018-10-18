from __future__ import print_function
import IMP
import IMP.container
import IMP.core
import IMP.display
import IMP.npctransport
import IMP.rmf
import IMP.test
import RMF
import math
import random

pradius=1
bottom = -5
top = 5
k_bounding= 10.0
k = 1.5 # in kal/mol/A
boxw= top * 4

class ZBiasTests(IMP.test.TestCase):
    def test_linear_radial_score(self):
        """Test radial force signleton score"""
        # setup model
        m= IMP.Model()
        p= IMP.Particle(m)
        d= IMP.core.XYZR.setup_particle(p)
        d.set_coordinates_are_optimized(True)
        d.set_radius(pradius)

        # Restraints list
        rs=[]

        # Bounding sphere restraint
        R=50
        bsphere= IMP.algebra.Sphere3D([0,0,0], R)
        bsss = IMP.core.BoundingSphere3DSingletonScore(IMP.core.HarmonicUpperBound(0,
                                                                                   k_bounding),
                                                       bsphere)
        rs.append(IMP.core.SingletonRestraint(m, bsss, p.get_index(), "bounding sphere"))

        # Radial restraint
        lrss= IMP.npctransport.LinearRadialSingletonScore(bsphere.get_center(),
                                                          k)
        lrsr= IMP.core.SingletonRestraint(m, lrss, p.get_index(), "radial")
        rs.append(lrsr)

        # Scoring function:
        sf= IMP.core.RestraintsScoringFunction(rs)

        # steep descent out
        cg= IMP.core.SteepestDescent(m)
        cg.set_scoring_function(sf)
        #        cg.set_log_level(IMP.VERBOSE)

        # RMF file optimizer state
        rmf = RMF.create_rmf_file("test_radial_ss.rmf")#output 1-1
        rmf.set_description("Gradient descent over radial singleton score test")
        h= IMP.atom.Hierarchy.setup_particle(p)
        IMP.rmf.add_hierarchy(rmf, h)
        IMP.rmf.add_restraints(rmf, rs)
        IMP.rmf.add_geometry(rmf, IMP.display.SphereGeometry(bsphere))
        sos = IMP.rmf.SaveOptimizerState(m, rmf)
        sos.update_always("initial conformation")
        sos.set_log_level(IMP.SILENT)
#        sos.set_simulator(cg)
        cg.add_optimizer_state(sos)
        print("RMF file is", rmf.get_name())

        # position d randonmly within excluded zrange
        d.set_coordinates([0.1,0.1,0.1])
        print("BEFORE=", d.get_coordinates())
        for i in range(0,5000):
            s=cg.optimize(200)
            cur_R= (d.get_coordinates()-bsphere.get_center()).get_magnitude()
            print("{:7d}".format(i), d.get_coordinates(), cur_R, s)
            self.assertAlmostEqual( s, - k * cur_R, delta=0.5 )

        print("AFTER=", d.get_coordinates(),
              "R =", cur_R,
              "Score =", s)
        self.assert_( cur_R > bsphere.get_radius() - pradius - 0.2 )


if __name__ == '__main__':
    IMP.test.main()
