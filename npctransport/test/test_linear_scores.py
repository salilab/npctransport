import IMP
import IMP.test
import IMP.npctransport
import IMP.container
import math

radius=3

class Tests(IMP.test.TestCase):
    def _create_particle(self, m):
        p= IMP.Particle(m)
        d= IMP.core.XYZR.setup_particle(p)
        d.set_radius(radius)
        d.set_coordinates_are_optimized(True)
        return d
    def test_repulsion(self):
        """Check linear soft sphere scores"""
        m= IMP.Model()
        m.set_log_level(IMP.SILENT)
        ds= [self._create_particle(m) for i in range(0,2)]
        ps= IMP.npctransport.LinearSoftSpherePairScore(1)
        ds[0].set_coordinates(IMP.algebra.Vector3D(0,0,0));
        ds[1].set_coordinates(IMP.algebra.Vector3D(0,0,0))
        ps.set_log_level(IMP.SILENT)
        sf= IMP.core.PairRestraint(ps, ds).create_scoring_function();
        for i in range(1,100):
            ds[1].set_coordinate(0, .1*i)
            e= sf.evaluate(True)
            deriv=ds[1].get_derivatives()
            print .1*i, e, deriv
            if .1*i > 2*radius:
                self.assertEqual(e, 0)
            elif .1*i < 2*radius:
                self.assert_(deriv[0] < 0)

    def test_interaction(self):
        """Check linear soft sphere scores"""
        m= IMP.Model()
        m.set_log_level(IMP.SILENT)
        ds= [self._create_particle(m) for i in range(0,2)]
        ps= IMP.npctransport.LinearInteractionPairScore(1, 2, .5)
        fps= IMP.npctransport.FunctorLinearInteractionPairScore(1, 2, .5)
        ds[0].set_coordinates(IMP.algebra.Vector3D(0,0,0));
        ds[1].set_coordinates(IMP.algebra.Vector3D(0,0,0))
        ps.set_log_level(IMP.SILENT)
        for i in range(0,100):
            ds[1].set_coordinate(0, .1*i)
            e= ps.evaluate(ds, None)
            ef=e= fps.evaluate(ds, None)
            print .1*i, e, ef

if __name__ == '__main__':
    IMP.test.main()
