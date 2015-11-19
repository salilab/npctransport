from __future__ import print_function
import IMP
import IMP.test
import IMP.npctransport
import math
import IMP.display
import random

radius=1
bottom = -5
top = 5
k = 10.0
boxw= top * 4

class ZBiasTests(IMP.test.TestCase):
    def test_z_bias_(self):
        """Check z-axis bias singleton score"""
        # setup model
        m= IMP.Model()
        p= IMP.Particle(m)
        d= IMP.core.XYZR.setup_particle(p)
        d.set_coordinates_are_optimized(True)
        d.set_radius(radius)
        bb= IMP.algebra.BoundingBox3D(IMP.algebra.Vector3D(-boxw, -boxw, -boxw),
                                      IMP.algebra.Vector3D(boxw, boxw, boxw))
        zbss= \
            IMP.npctransport.ZBiasSingletonScore(k)
        sr= IMP.core.SingletonRestraint(m, zbss, p.get_index(), "zbias")

        # position d randonmly within excluded zrange
        d.set_coordinates([1,1,100])
        print("BEFORE=", d.get_coordinates())
        # steep descent out of excluded zone (hopefully)
        cg= IMP.core.SteepestDescent(m)
        cg.set_scoring_function(sr)
        cg.set_log_level(IMP.VERBOSE)
        for i in range(0,10000):
            s=cg.optimize(3)
#            print d.get_coordinates(), s
            if s <= 0:
                break
        print("AFTER=", d.get_coordinates(), "Score =", s)
        self.assert_(d.get_z() <= 0.01 and
                     d.get_x() == 1 and
                     d.get_y() == 1)
if __name__ == '__main__':
    IMP.test.main()
