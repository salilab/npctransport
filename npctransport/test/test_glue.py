import IMP
import IMP.test
import IMP.npctransport
import RMF
import IMP.container
import math

radius=5

class Tests(IMP.test.TestCase):
    def _create_particle(self, m):
        p= IMP.Particle(m)
        d= IMP.core.XYZR.setup_particle(p)
        d.set_radius(radius)
        IMP.atom.Hierarchy.setup_particle(p)
        IMP.atom.Mass.setup_particle(p, 1.0)
        IMP.core.RigidBody.setup_particle(p, IMP.algebra.ReferenceFrame3D())
        d.set_coordinates_are_optimized(True)
        IMP.atom.RigidBodyDiffusion.setup_particle(p)
        return d
    def _test_one(self, site_range, site_k, nonspec_range, nonspec_k,
                  soft_sphere_k, dt):
        m= IMP.Model()
        m.set_log_level(IMP.SILENT)
        ds= [self._create_particle(m) for i in range(0,2)]
        ds[0].set_coordinates(IMP.algebra.Vector3D(0,0,0))
        ds[1].set_coordinates(IMP.algebra.Vector3D(2*radius,0,0))
        types=[IMP.core.ParticleType(d.get_name()+" type") for d in ds]
        for d in zip(types, ds):
          IMP.core.Typed.setup_particle(d[1], d[0])
        sites=([IMP.algebra.Vector3D(radius, 0,0)], [IMP.algebra.Vector3D(-radius, 0,0)])
        ps= IMP.npctransport.SitesPairScore(site_range, site_k, nonspec_range,
                                            nonspec_k, soft_sphere_k, sites[0], sites[1])
        ps.set_log_level(IMP.VERBOSE)
        r= IMP.core.PairRestraint(ps, ds)
        m.add_restraint(r)
        bd= IMP.atom.BrownianDynamics(m)
        bd.set_maximum_time_step(dt)
        f= RMF.create_rmf_file(self.get_tmp_file_name("glue.rmf"))
        for d in zip(types, sites):
          IMP.npctransport.add_sites(f, d[0], .5*radius, d[1])
        w= IMP.npctransport.add_hierarchies(f, ds)
        sos= IMP.rmf.SaveOptimizerState(f)
        bd.add_optimizer_state(sos)
        bd.optimize(1000)
    def test_cone_construction2(self):
        """Check interaction score repulsion"""
        k=10
        dt=IMP.npctransport.get_time_step(1, k, radius)
        self._test_one(radius, k, .2*radius, .5*k, k, dt)


if __name__ == '__main__':
    IMP.test.main()
