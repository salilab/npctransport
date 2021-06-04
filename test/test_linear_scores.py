from __future__ import print_function
import IMP
import IMP.test
import IMP.npctransport
import IMP.container
import IMP.score_functor
import math

radius=3

class Tests(IMP.test.TestCase):
    def _create_xyzr(self, m):
        p= IMP.Particle(m)
        d= IMP.core.XYZR.setup_particle(p)
        d.set_radius(radius)
        d.set_coordinates_are_optimized(True)
        return d
    def test_repulsion(self):
        """Check linear soft sphere scores"""
        m= IMP.Model()
        m.set_log_level(IMP.SILENT)
        ds= [self._create_xyzr(m) for i in range(0,2)]
        dsi=[x.get_particle_index() for x in ds]
        ps= IMP.npctransport.LinearSoftSpherePairScore(1)
        ds[0].set_coordinates(IMP.algebra.Vector3D(0,0,0));
        ds[1].set_coordinates(IMP.algebra.Vector3D(0,0,0))
        ps.set_log_level(IMP.SILENT)
        sf= IMP.core.PairRestraint(m, ps, dsi).create_scoring_function();
        for i in range(1,100):
            ds[1].set_coordinate(0, .1*i)
            e= sf.evaluate(True)
            deriv=ds[1].get_derivatives()
            print(.1*i, e, deriv)
            if .1*i > 2*radius:
                self.assertEqual(e, 0)
            elif .1*i < 2*radius:
                self.assertLess(deriv[0], 0)

    def test_interaction(self):
        """Check linear interaction scores"""
        m= IMP.Model()
        m.set_log_level(IMP.SILENT)
        ds= [self._create_xyzr(m) for i in range(0,2)]
        dsi=[x.get_particle_index() for x in ds]
        ps= IMP.npctransport.LinearInteractionPairScore(1, 2, .5)
        fps= IMP.npctransport.FunctorLinearInteractionPairScore(1, 2, .5)
        ds[0].set_coordinates(IMP.algebra.Vector3D(0,0,0));
        ds[1].set_coordinates(IMP.algebra.Vector3D(0,0,0))
        ps.set_log_level(IMP.SILENT)
        for i in range(0,100):
            ds[1].set_coordinate(0, .1*i)
#            IMP.score_functor.evaluate_index(ps, m, dsi, None)
            e= ps.evaluate_index(m,dsi, None)
            ef=e= fps.evaluate_index(m,dsi, None)
            print(.1*i, e, ef)

if __name__ == '__main__':
    IMP.test.main()
