import IMP
import IMP.test
import IMP.npctransport
import IMP.container
import IMP.rmf
import RMF
import math

radius=5

class ConeTests(IMP.test.TestCase):
    def _create_particle(self, m):
        p= IMP.Particle(m)
        d= IMP.core.XYZR.setup_particle(p)
        d.set_radius(radius)
        d.set_coordinates_are_optimized(True)
        h= IMP.atom.Hierarchy.setup_particle(p)
        m= IMP.atom.Mass.setup_particle(p, 1)
        diff= IMP.atom.Diffusion.setup_particle(p)
        return d
    def _randomize(self, ds, bb):
        for d in ds:
            d.set_coordinates(IMP.algebra.get_random_vector_in(bb))
    def _show(self, ds, w):
        for d in ds:
            g= IMP.core.XYZRGeometry(d);
            w.add_geometry(g)
    def test_cone_construction(self):
        """Check linear well"""
        m= IMP.Model()
        m.set_log_level(IMP.SILENT)
        ds= [self._create_particle(m) for i in range(0,2)]
        ds[1].set_coordinates(IMP.algebra.Vector3D(0,2*radius,0))
        ps= IMP.npctransport.LinearWellPairScore(2*radius, 20)
        r= IMP.core.PairRestraint(ps, ds)
        bd= IMP.atom.BrownianDynamics(m)
        bd.set_maximum_time_step(100)
        bd.set_scoring_function([r])
        f= RMF.create_rmf_file(self.get_tmp_file_name("well.rmf"))
        IMP.rmf.add_hierarchies(f, ds)
        IMP.rmf.add_restraints(f, [r])
        w= IMP.rmf.SaveOptimizerState(f)
        bd.add_optimizer_state(w)
        bd.optimize(1000)
        dist= IMP.core.get_distance(ds[0], ds[1])
        print dist
        self.assertAlmostEqual(dist, 0, delta=radius)

if __name__ == '__main__':
    IMP.test.main()
