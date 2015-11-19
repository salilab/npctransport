from __future__ import print_function
import IMP
import IMP.test
import IMP.npctransport
import math
import IMP.display
import random

radius=1
#random.uniform(1,12)
slab_radius=3
#random.uniform(radius+2, radius+15)
slab_height=5
#random.uniform(5,30)
boxw= max([1.5*slab_radius,slab_height])
def out_slab(d):
    c= d.get_coordinates()
    if c[2]> slab_height/2.0+radius-.5:
        return True
    if c[2]< -slab_height/2.0-radius+.5:
        return True
    rxz= (c[0]**2+c[1]**2)**.5
    if rxz +radius < slab_radius+.5:
        return True
    print(rxz, slab_radius-radius, c[2], slab_height/2.0+radius)
    return False
class ConeTests(IMP.test.TestCase):
    def test_cone_construction(self):
        """Check slab singleton score"""
        print("radius", radius, "slab radius", slab_radius, "slab_height", slab_height)
        m= IMP.Model()
        p= IMP.Particle(m)
        d= IMP.core.XYZR.setup_particle(p)
        d.set_coordinates_are_optimized(True)
        d.set_radius(radius)
        bb= IMP.algebra.BoundingBox3D(IMP.algebra.Vector3D(-boxw, -boxw, -boxw),
                                      IMP.algebra.Vector3D(boxw,boxw,boxw))
        slabss= IMP.npctransport.SlabSingletonScore(slab_height, slab_radius, 1)
        r= IMP.core.SingletonRestraint(m, slabss, p.get_index(), "slab")
        while out_slab(d):
            d.set_coordinates(IMP.algebra.get_random_vector_in(bb))
        print(d.get_coordinates())
        w= IMP.display.PymolWriter(self.get_tmp_file_name("slabss.pym"))
        w.set_frame(0)
        g=IMP.core.XYZRGeometry(d)
        sg= IMP.npctransport.SlabWireGeometry(slab_height, slab_radius, boxw)
        w.add_geometry([g, sg])
        cg= IMP.core.SteepestDescent(m)
        cg.set_scoring_function(r)
        cg.set_log_level(IMP.VERBOSE)
        for i in range(0,100):
            s=cg.optimize(10)
            w.set_frame(i+1)
            w.add_geometry([g, sg])
            if s==0:
                break
        print(d.get_coordinates())
        self.assert_(out_slab(d))
if __name__ == '__main__':
    IMP.test.main()
