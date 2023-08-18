from __future__ import print_function
import IMP
import IMP.test
import IMP.npctransport
import math
import IMP.display
import random

debug=True
radius=1
slab_pore_radius=5
slab_height=3
boxw= 2*max([3*slab_pore_radius,slab_height])
def out_slab(p, slab):
    ''' verify particle p is out of slab.
        p is assumed to be decorated by XYZR '''
    d=IMP.core.XYZR(p)
    c= d.get_coordinates()
    if c[2]> slab_height/2.0+radius-.2:
        return True
    if c[2]< -slab_height/2.0-radius+.2:
        return True
    rxy= (c[0]**2+c[1]**2)**.5
    print("out_slab() - c: ", c, " rxy,max_rxy: ", rxy,
          slab.get_pore_radius()-radius, " z,min_z:", c[2], slab_height/2.0+radius)
    print("rxy+radius", rxy+radius)
    print("pore radius + .2 = ", slab.get_pore_radius()+.2)
    IS_OUT = (rxy + radius) < (slab.get_pore_radius() + .2)
    print("(rxy + radius) < (slab.get_pore_radius() + .2) = ", IS_OUT)
    if IS_OUT:
        print("out_slab returning true")
        return True
    print("out_slab is False - c: ", c, " rxy,max_rxy: ", rxy,
          slab.get_pore_radius()-radius, " z,min_z:", c[2], slab_height/2.0+radius)
    return False


class CylindricalPoreSSTest(IMP.test.TestCase):
    def _test_optimization(self, pymol_fname=None):
        '''
        pymol_fname - file name for output (stored in tmp if debug==True)
        '''
        global debug
        d= self.d # diffusing particle
        slab= self.slab
        opt= self.opt # optimizer
        if not debug:
            pymol_fname= self.get_tmp_file_name(pymol_fname)
        # Make pymol writer and geometries:
        if pymol_fname is not None:
            w= IMP.display.PymolWriter(pymol_fname)
            w.set_frame(0)
            sg= IMP.npctransport.SlabWithCylindricalPoreWireGeometry(slab_height, slab_pore_radius, boxw)
            g=IMP.core.XYZRGeometry(d)
            w.add_geometry([g, sg])
        else:
            w= None
        # Optimize:
        if(debug): print(d.get_coordinates())
        for i in range(0,5000):
            s=opt.optimize(1)
            if w is not None:
                w.set_frame(i+1)
                w.add_geometry([g, sg])
            if(debug):
                print("Score=%.4f at i=%d" % (s,i), " d=", d.get_coordinates(),
                      " pore radius=",slab.get_pore_radius())
                print("Derivative XYZ", IMP.core.XYZ(d).get_derivatives())
                print("Pore Radius derivative",
                      slab.get_particle().get_derivative
                      ( IMP.npctransport.SlabWithPore.get_pore_radius_key() ) )
            if abs(s-0)<0.0001:
                if(debug):
                    print("*** BREAKING ***")
                break
        print("Final coordinates: ", d.get_coordinates(),  " score ", s, "Pore radius", slab.get_pore_radius())
        OUT_SLAB = out_slab(d,slab)
        print("OUT_SLAB = ", OUT_SLAB)
        self.assertTrue(OUT_SLAB)

    def _initialize_model(self):
        print("radius", radius, "slab radius", slab_pore_radius, "slab_height", slab_height)
        m=IMP.Model()
        self.m= m
        p= IMP.Particle(m,"diffuser")
        d= IMP.core.XYZR.setup_particle(p)
        self.d= d
        d.set_coordinates_are_optimized(True)
        d.set_radius(radius)
        bb= IMP.algebra.BoundingBox3D(0.5*IMP.algebra.Vector3D(-boxw, -boxw, -boxw),
                                      0.5*IMP.algebra.Vector3D(boxw,boxw,boxw))
        p_slab= IMP.Particle(m, "slab")
        IMP.npctransport.SlabWithCylindricalPore.setup_particle \
            (p_slab, slab_height, slab_pore_radius)
        self.assertTrue(
            IMP.npctransport.SlabWithCylindricalPore.get_is_setup(p_slab))
        # test cast to slab
        slab= IMP.npctransport.SlabWithPore(p_slab)
        self.slab= slab
        self.assertEqual(slab.get_pore_radius(),slab_pore_radius)
        self.assertEqual(slab.get_thickness(),slab_height);
        slabps= IMP.npctransport.SlabWithCylindricalPorePairScore(1.0)
        self.assertFalse(slab.get_pore_radius_is_optimized()) # verify correct default value
        r= IMP.core.PairRestraint(m, slabps, [p_slab.get_index(), p.get_index()], "slab")
        self.r= r
        score_init= slabps.evaluate_index(m,
                                          [p_slab.get_index(), p.get_index()],
                                          IMP.DerivativeAccumulator(1.0))
        self.assertEqual(score_init,0.0)
        # Create optimizer:
        cg= IMP.core.SteepestDescent(m)
        cg.set_scoring_function(r)
        cg.set_step_size(0.01)
        self.opt= cg

    def test_slab_pair_score(self):
        """Check slab pair score"""
        self._initialize_model()
        d= self.d # diffusing particle
        slab= self.slab

        print("\n== (I) Testing non-optimizable pore radius ==")

        print("\nTest from random position")
        bb_half= IMP.algebra.BoundingBox3D(0.25*IMP.algebra.Vector3D(-boxw, -boxw, -boxw),
                                           0.25*IMP.algebra.Vector3D(boxw,boxw,boxw))
        while out_slab(d,slab):
            d.set_coordinates(IMP.algebra.get_random_vector_in(bb_half))
        self._test_optimization('tmp.pym')
        self.assertEqual(slab_pore_radius, slab.get_pore_radius())

        print("\nTest from vertical position")
        d.set_coordinates([slab_pore_radius+2*radius,0,0.4])
        self.assertFalse(out_slab(d,slab))
        self._test_optimization('tmp3.pym')
        self.assertEqual(slab_pore_radius, slab.get_pore_radius())

        print("\nTest from horizontal position")
        d.set_coordinates([slab_pore_radius+radius,0,0])
        self.assertFalse(out_slab(d,slab))
        self._test_optimization('tmp2.pym')
        self.assertEqual(slab_pore_radius, slab.get_pore_radius())

        print("\n== (II) Testing optimizable pore radius ==")

        slab.set_pore_radius_is_optimized(True)
        print("\nTest pore radius from vertical position")
        d.set_coordinates([slab_pore_radius+2*radius,0,0.4])
        self._test_optimization('tmp4.pym')
        self.assertEqual(slab_pore_radius, slab.get_pore_radius()) # vertical position shouldn't have affected pore radius

        print("\nTest pore radius from horizontal position")
        d.set_coordinates([slab_pore_radius+radius,0,0.1])
        self._test_optimization('tmp5.pym')
        self.assertAlmostEqual(slab.get_pore_radius(), 6.01, delta=0.025)

    def test_pore_radius_score_with_slab_pair_score(self):
        '''Test combination of pore radius and slab scores'''
        self._initialize_model()
        m= self.m
        d=  self.d
        slab= self.slab

        # append pore radius score to existing restraint
        prss=IMP.npctransport.PoreRadiusSingletonScore(slab_pore_radius, 1.0)
        r_prss= IMP.core.SingletonRestraint(m,
                                            prss,
                                            slab.get_particle().get_index(),
                                            "pore radius score")
        self.opt.set_scoring_function(IMP.core.RestraintsScoringFunction([self.r, r_prss]))
        slab.set_pore_radius_is_optimized(True)
        # optimize with new score
        print("Testing with pore radius k= 1.0")
        d.set_coordinates([slab_pore_radius+radius,0,0.1])
        self._test_optimization('tmp5.pym')
        self.assertAlmostEqual(slab.get_pore_radius(), 5.37, delta=0.025)
        # optimize with larger k
        print("Testing with pore radius k= 10.0")
        prss.set_k (10.0)
        d.set_coordinates([slab_pore_radius+radius,0,0.1])
        self._test_optimization('tmp5.pym')
        self.assertAlmostEqual(slab.get_pore_radius(), 5.04, delta=0.025)
        # optimize with tiny k
        print("Testing with pore radius k= 0.001")
        prss.set_k (0.001)
        d.set_coordinates([slab_pore_radius+radius,0,0.1])
        self._test_optimization('tmp5.pym')
        self.assertAlmostEqual(slab.get_pore_radius(), 6.02, delta=0.025)

if __name__ == '__main__':
    IMP.test.main()
