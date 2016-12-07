from __future__ import print_function
import IMP
import IMP.test
import IMP.npctransport
import math
import IMP.display
import random

radius=1
#random.uniform(1,12)
slab_radius=5
#random.uniform(radius+2, radius+15)
slab_height=3
#random.uniform(5,30)
boxw= 2*max([3*slab_radius,slab_height])
def out_slab(p):
    ''' verify particle p is out of slab.
        p is assumed to be decorated by XYZR '''
    d=IMP.core.XYZR(p)
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
    def test_slab_pair_score(self):
        """Check slab pair score"""
        print("radius", radius, "slab radius", slab_radius, "slab_height", slab_height)
        m= IMP.Model()
        p= IMP.Particle(m,"diffuser")
        d= IMP.core.XYZR.setup_particle(p)
        d.set_coordinates_are_optimized(True)
        d.set_radius(radius)
        bb= IMP.algebra.BoundingBox3D(0.5*IMP.algebra.Vector3D(-boxw, -boxw, -boxw),
                                      0.5*IMP.algebra.Vector3D(boxw,boxw,boxw))
        bb_half= IMP.algebra.BoundingBox3D(0.25*IMP.algebra.Vector3D(-boxw, -boxw, -boxw),
                                      0.25*IMP.algebra.Vector3D(boxw,boxw,boxw))
        p_slab= IMP.Particle(m, "slab")
        slab_orig= IMP.npctransport.SlabWithCylindricalPore.setup_particle \
                   (p_slab, slab_height, slab_radius)
        self.assert_(IMP.npctransport.SlabWithPore.get_is_setup(p_slab))
        # test cast to slab
        slab= IMP.npctransport.SlabWithPore(p_slab)
        self.assertEqual(slab.get_pore_radius(),slab_radius)
        self.assertEqual(slab.get_thickness(),slab_height);
        slabps= IMP.npctransport.SlabWithCylindricalPorePairScore(1.0)
        r= IMP.core.PairRestraint(m, slabps, [p_slab.get_index(), p.get_index()], "slab")
        score_init= slabps.evaluate_index(m,
                                          [p_slab.get_index(), p.get_index()],
                                          IMP.DerivativeAccumulator(1.0))
        self.assertEqual(score_init,0.0)
        while out_slab(p):
            d.set_coordinates(IMP.algebra.get_random_vector_in(bb_half))
        print(d.get_coordinates())
        w= IMP.display.PymolWriter("tmp.pym") #self.get_tmp_file_name("slabps.pym"))
        w.set_frame(0)
        g=IMP.core.XYZRGeometry(d)
        sg= IMP.npctransport.SlabWithCylindricalPoreWireGeometry(slab_height, slab_radius, boxw)
        w.add_geometry([g, sg])
        cg= IMP.core.SteepestDescent(m)
        cg.set_scoring_function(r)
        cg.set_log_level(IMP.VERBOSE)
        for i in range(0,1000):
            s=cg.optimize(1)
            w.set_frame(i+1)
            w.add_geometry([g, sg])
            if s==0:
                break
        print(d.get_coordinates())
        self.assert_(out_slab(d))
if __name__ == '__main__':
    IMP.test.main()
