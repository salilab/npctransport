from __future__ import print_function
import IMP
import IMP.test
import IMP.algebra
import IMP.core
import IMP.npctransport
import math
import IMP.display
import random

radius=1
#random.uniform(1,12)
#random.uniform(radius+2, radius+15)
slab_thickness=5
slab_radius=5+0.5*slab_thickness
#random.uniform(5,30)
boxw= 2*max([1*slab_radius,slab_thickness])

def out_slab(p, ALLOWED_OVERLAP=0.0):
    '''
    verify that spherical (XYZR) particle p is out of slab with
    thickness slab_thickness, and toroidal pore with major radius
    slab_radius and minor radius 0.5*slab_thickness.

    penetration of ALLOWED_OVERLAP is allowed (would still result in True). Use negative value to force extra slack
    '''
    EPS=0.000001
    R=slab_radius # major radius of torus
    r=0.5*slab_thickness # minor radius of torus
    d= IMP.core.XYZR(p)
    print("out_slab begins: ", d,R,r)
    xyz= d.get_coordinates()
    if xyz[2] - d.get_radius() + ALLOWED_OVERLAP > r:
        print("Above slab")
        return True
    if xyz[2] + d.get_radius() - ALLOWED_OVERLAP < -r:
        print("Below slab")
        return True
    print("In slab vertically")
    xy0= IMP.algebra.Vector3D(xyz[0], xyz[1], 0.0)
    rxy0= xy0.get_magnitude()
    print("xy0", xy0, "magnitude", rxy0)
    if rxy0 + d.get_radius() > R:
        print("Out of outer radius perimeter (for entire sphere)")
        return False
    print("Inside outer radius perimeter")
    if(rxy0>EPS):
        xy0_major= xy0 * ( R / rxy0 ) # same direction as xy0 but on major radius
    else:
        xy0_major=IMP.algebra.Vector3D(0, R, 0.0)
    print("xy0_major", xy0_major)
    distance= (xyz - xy0_major).get_magnitude()
    print("Distance", distance)
    is_overlap = (distance - d.get_radius() <= r - ALLOWED_OVERLAP)
    print("Is_overlap", is_overlap)
    return not is_overlap

class ConeTests(IMP.test.TestCase):
    def test_slab_singleton_score(self):
        """Check slab singleton score"""
        IMP.set_log_level(IMP.SILENT)
        print("radius", radius, "slab radius", slab_radius, "slab_thickness", slab_thickness)
        m= IMP.Model()
        m.set_log_level(IMP.SILENT)
        p= IMP.Particle(m)
        d= IMP.core.XYZR.setup_particle(p)
        d.set_coordinates_are_optimized(True)
        d.set_radius(radius)
        bb= IMP.algebra.BoundingBox3D(0.5*IMP.algebra.Vector3D(-boxw, -boxw, -boxw),
                                      0.5*IMP.algebra.Vector3D(boxw,boxw,boxw))
        slabss= IMP.npctransport.SlabWithToroidalPoreSingletonScore \
                (slab_thickness, slab_radius, 1)
        slabss.set_log_level(IMP.SILENT)
        self.assertEqual(slabss.get_bottom_z(),-0.5*slab_thickness);
        self.assertEqual(slabss.get_top_z(),+0.5*slab_thickness);
        r= IMP.core.SingletonRestraint(m, slabss, p.get_index(), "slab")
        while out_slab(p, ALLOWED_OVERLAP=-0.5):
            d.set_coordinates(IMP.algebra.get_random_vector_in(bb))
        print("Penetrating slab: ", d.get_coordinates())
        pym_fname="tmp.pym" #self.get_tmp_file_name("slabss.pym")
        w= IMP.display.PymolWriter(pym_fname)
        w.set_frame(0)
        g=IMP.core.XYZRGeometry(d)
        sg= IMP.npctransport.SlabWithToroidalPoreWireGeometry(slab_thickness, slab_radius, boxw)
        w.add_geometry([g, sg])
        cg= IMP.core.SteepestDescent(m)
        cg.set_scoring_function(r)
        cg.set_log_level(IMP.SILENT)
        f=1
        for i in range(0,1000):
            s=cg.optimize(1)
            if(i % 10 == 0):
                w.set_frame(f)
                f=f+1
                w.add_geometry([g, sg])
            print(i, d.get_coordinates())
            if s==0:
                break
        print(d.get_coordinates())
        self.assert_(out_slab(d, ALLOWED_OVERLAP=0.05))
if __name__ == '__main__':
    IMP.test.main()
