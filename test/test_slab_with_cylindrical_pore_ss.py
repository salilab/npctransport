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
    if c[2]> slab_height/2.0+radius-.1:
        return True
    if c[2]< -slab_height/2.0-radius+.1:
        return True
    rxy= (c[0]**2+c[1]**2)**.5
    if rxy +radius < slab_radius+.1:
        return True
    print("out_slab is False - c: ", c, " rxy,max_rxy: ", rxy, slab_radius-radius, " z,min_z:", c[2], slab_height/2.0+radius)
    return False


class CylindricalPoreSSTest(IMP.test.TestCase):
    def _test_optimization(self, m, d, opt, pymol_fname=None):
        '''
        m - model
        d - diffusing particle
        opt - optimizer
        '''
        debug=False
        if not debug:
            pymol_fname= self.get_tmp_file_name(pymol_fname)
        # Make pymol writer and geometries:
        if pymol_fname is not None:
            w= IMP.display.PymolWriter(pymol_fname)
            w.set_frame(0)
            sg= IMP.npctransport.SlabWithCylindricalPoreWireGeometry(slab_height, slab_radius, boxw)
            g=IMP.core.XYZRGeometry(d)
            w.add_geometry([g, sg])
        else:
            w= None
        # Optimize:
        if(debug): print(d.get_coordinates())
        for i in range(0,2000):
            s=opt.optimize(1)
            if w is not None:
                w.set_frame(i+1)
                w.add_geometry([g, sg])
            if(debug): print("Score=%.2f at i=%d" % (s,i), d.get_coordinates())
            if s==0:
                break
        print("Final coordinates: ", d.get_coordinates(),  " score ", s)
        self.assert_(out_slab(d))


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
        IMP.npctransport.SlabWithCylindricalPore.setup_particle \
            (p_slab, slab_height, slab_radius)
        self.assert_(IMP.npctransport.SlabWithCylindricalPore.get_is_setup(p_slab))
        # test cast to slab
        slab= IMP.npctransport.SlabWithPore(p_slab)
        self.assertEqual(slab.get_pore_radius(),slab_radius)
        self.assertEqual(slab.get_thickness(),slab_height);
        slabps= IMP.npctransport.SlabWithCylindricalPorePairScore(1.0)
        self.assert_(not IMP.npctransport.SlabWithPore(p_slab).get_pore_radius_is_optimized()) # verify correct default value
        r= IMP.core.PairRestraint(m, slabps, [p_slab.get_index(), p.get_index()], "slab")
        score_init= slabps.evaluate_index(m,
                                          [p_slab.get_index(), p.get_index()],
                                          IMP.DerivativeAccumulator(1.0))
        self.assertEqual(score_init,0.0)

        # Create optimizer:
        cg= IMP.core.SteepestDescent(m)
        cg.set_scoring_function(r)
        cg.set_log_level(IMP.VERBOSE)
        IMP.add_to_log("hello\n")

        slabps.set_log_level(IMP.VERBOSE)
        IMP.set_log_level(IMP.VERBOSE)
        r.set_log_level(IMP.VERBOSE)


        print("Test from random position")
        while out_slab(p):
            d.set_coordinates(IMP.algebra.get_random_vector_in(bb_half))
        self._test_optimization(m, d, cg, 'tmp.pym')

        print("Test from horizontal position")
        d.set_coordinates([slab_radius+radius,0,0])
        self.assert_(not out_slab(d))
        self._test_optimization(m, d, cg, 'tmp2.pym')

        print("Test from vertical position")
        d.set_coordinates([slab_radius+2*radius,0,0.4])
        self.assert_(not out_slab(d))
        self._test_optimization(m, d, cg, 'tmp3.pym')

if __name__ == '__main__':
    IMP.test.main()
