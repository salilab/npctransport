from __future__ import print_function
import IMP
import IMP.test
import IMP.npctransport
import IMP.container
import math

radius=5

class ConeTests(IMP.test.TestCase):
    def _create_particle(self, m):
        p= IMP.Particle(m)
        d= IMP.core.XYZR.setup_particle(p)
        d.set_radius(radius)
        d.set_coordinates_are_optimized(True)
        return d
    def _randomize(self, ds, bb):
        for d in ds:
            d.set_coordinates(IMP.algebra.get_random_vector_in(bb))
    def _show(self, ds, w):
        for d in ds:
            g= IMP.core.XYZRGeometry(d);
            w.add_geometry(g)
    def test_cone_construction(self):
        """Check linear soft sphere"""
        m= IMP.Model()
        m.set_log_level(IMP.SILENT)
        ds= [self._create_particle(m) for i in range(0,2)]
        apps= IMP.container.AllPairContainer(IMP.container.ListSingletonContainer(m, IMP.get_indexes(ds)))
        ps= IMP.npctransport.LinearSoftSpherePairScore(10)
        ps.set_log_level(IMP.VERBOSE)
        m.set_log_level(IMP.SILENT)
        r= IMP.container.PairsRestraint(ps, apps)
        sf = IMP.core.RestraintsScoringFunction([r])
        bb= IMP.algebra.get_cube_3d(5)
        ds[0].set_coordinates(IMP.algebra.Vector3D(0,0,0))
        ds[0].set_coordinates_are_optimized(False)
        self._randomize(ds[1:], bb)
        w= IMP.display.PymolWriter(self.get_tmp_file_name("linearss.pym"))
        w.set_frame(0)
        self._show(ds, w)
        cg= IMP.core.SteepestDescent(m)
        cg.set_scoring_function(sf)
        cg.set_log_level(IMP.VERBOSE)
        cg.set_step_size(.01)
        cg.set_maximum_step_size(.01)
        print("initial", sf.evaluate(True))
        cg.optimize(1000)
        w.set_frame(1)
        self._show(ds, w)
        for d in ds:
            for d1 in ds:
                if d!=d1:
                    print(d, d1)
                    self.assert_(IMP.core.get_distance(d, d1) > -.1)
    def test_cone_construction2(self):
        """Check interaction score repulsion"""
        m= IMP.Model()
        rng=3
        m.set_log_level(IMP.SILENT)
        ds= [self._create_particle(m) for i in range(0,2)]
        ds[0].set_coordinates(IMP.algebra.Vector3D(0,0,0))
        ds[0].set_coordinates_are_optimized(False)
        bb= IMP.algebra.get_cube_3d(3)
        self._randomize(ds, bb)
        apps= IMP.container.AllPairContainer(IMP.container.ListSingletonContainer(m, IMP.get_indexes(ds)))
        ps= IMP.npctransport.LinearInteractionPairScore(radius*2, 10, rng)
        ps.set_log_level(IMP.VERBOSE)
        r= IMP.container.PairsRestraint(ps, apps)
        sf = IMP.core.RestraintsScoringFunction([r])
        bb= IMP.algebra.get_cube_3d(5)
        self._randomize(ds[1:], bb)
        w= IMP.display.PymolWriter(self.get_tmp_file_name("linearint.pym"))
        w.set_frame(0)
        self._show(ds, w)
        cg= IMP.core.SteepestDescent(m)
        cg.set_scoring_function(sf)
        cg.set_step_size(.003)
        cg.set_threshold(-10000)
        cg.set_maximum_step_size(.003)
        print("initial", sf.evaluate(True))
        cg.optimize(1000)
        w.set_frame(1)
        self._show(ds, w)
        for d in ds:
            for d1 in ds:
                if d!=d1:
                    print(d, d1)
                    self.assert_(IMP.core.get_distance(d, d1) > -.1)
        ss= IMP.algebra.Sphere3D(IMP.algebra.Vector3D(0,0,0),
                                                      radius+rng*.9)
        ds[1].set_coordinates(IMP.algebra.get_random_vector_on(ss))
        w.set_frame(2)
        self._show(ds, w)
        print("initial", sf.evaluate(True))
        cg.optimize(1000)
        w.set_frame(3)
        self._show(ds, w)
        for d in ds:
            for d1 in ds:
                if d!=d1:
                    print(d, d1)
                    self.assert_(IMP.core.get_distance(d, d1) > -.1\
                                   and IMP.core.get_distance(d, d1) < .1)

if __name__ == '__main__':
    IMP.test.main()
