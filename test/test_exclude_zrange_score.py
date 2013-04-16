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

def out_zrange(d):
    c= d.get_coordinates()
    return (c[2] > top or c[2] < bottom)
class ExcludeZRangeTests(IMP.test.TestCase):
    def test_cone_(self):
        """Check exclude z-range singleton score"""
        print "top", top, "bottom", bottom, "k", k, "particle-radius", radius
        # setup model
        m= IMP.Model()
        p= IMP.Particle(m)
        d= IMP.core.XYZR.setup_particle(p)
        d.set_coordinates_are_optimized(True)
        d.set_radius(radius)
        bb= IMP.algebra.BoundingBox3D(IMP.algebra.Vector3D(-boxw, -boxw, -boxw),
                                      IMP.algebra.Vector3D(boxw, boxw, boxw))
        bb_exclude = \
            IMP.algebra.BoundingBox3D(IMP.algebra.Vector3D(-boxw, -boxw, bottom),
                                      IMP.algebra.Vector3D(boxw, boxw, top))
        # create restrain on z-range exclusion
        exclude_zrange_ss= \
            IMP.npctransport.ExcludeZRangeSingletonScore(bottom, top, k)
        r= IMP.core.SingletonRestraint(exclude_zrange_ss, p, "slab")

        # position d randonmly within excluded zrange
        d.set_coordinates(IMP.algebra.get_random_vector_in(bb_exclude))
        print d.get_coordinates()
        # steep descent out of excluded zone (hopefully)
        cg= IMP.core.SteepestDescent(m)
        cg.set_scoring_function(r)
        cg.set_log_level(IMP.VERBOSE)
        for i in range(0,10000):
            s=cg.optimize(3)
            print d.get_coordinates(), s
            if s==0:
                break
        print d.get_coordinates()
        self.assert_(out_zrange(d))
if __name__ == '__main__':
    IMP.test.main()
