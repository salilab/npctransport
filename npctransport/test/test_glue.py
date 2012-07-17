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
    def test_one(self):
        """Check interaction score repulsion"""
        k=500
        dt=IMP.npctransport.get_time_step(1, k, radius)
        self._test_one(radius, k, .2*radius, .5*k, k, dt)

    def _test_two(self, site_range, site_k, nonspec_range, nonspec_k,
                  soft_sphere_k, dt):
        m= IMP.Model()
        m.set_log_level(IMP.SILENT)
        ds= [self._create_particle(m) for i in range(0,3)]
        ds[0].set_coordinates(IMP.algebra.Vector3D(0,0,0))
        ds[1].set_coordinates(IMP.algebra.Vector3D(2*radius,0,0))
        ds[2].set_coordinates(IMP.algebra.Vector3D(0,0,2*radius))
        types=[IMP.core.ParticleType(d.get_name()+" type") for d in ds]
        for d in zip(types, ds):
          IMP.core.Typed.setup_particle(d[1], d[0])
        sites=([IMP.algebra.Vector3D(radius, 0,0),
                IMP.algebra.Vector3D(0, 0,radius)], [IMP.algebra.Vector3D(-radius, 0,0)],
               [IMP.algebra.Vector3D(0, 0,-radius)])
        rs=[]
        for p in [(0,1), (1,2), (0,2)]:
          ps= IMP.npctransport.SitesPairScore(site_range, site_k, nonspec_range,
                                              nonspec_k, soft_sphere_k, sites[p[0]], sites[p[1]])
          ps.set_log_level(IMP.VERBOSE)
          r=IMP.core.PairRestraint(ps, (ds[p[0]], ds[p[1]]))
          rs.append(r)
        bd= IMP.atom.BrownianDynamics(m)
        bd.set_scoring_function(rs)
        bd.set_maximum_time_step(dt)
        f= RMF.create_rmf_file(self.get_tmp_file_name("glue2.rmf"))
        for d in zip(types, sites):
          IMP.npctransport.add_sites(f, d[0], .5*radius, d[1])
        w= IMP.npctransport.add_hierarchies(f, ds)
        IMP.rmf.add_restraints(f, rs)
        sos= IMP.rmf.SaveOptimizerState(f)
        bd.add_optimizer_state(sos)
        bd.optimize(2000)
    def test_two(self):
        """Check two interactions"""
        k=500
        dt=IMP.npctransport.get_time_step(1, k, radius)
        self._test_two(radius, k, .2*radius, .5*k, k, dt)

    def _create_restraint_three(self, m, ds, site_range, site_k, nonspec_range, nonspec_k,
                                soft_sphere_k, f=None):
        sites=([IMP.algebra.Vector3D(radius, 0,0), IMP.algebra.Vector3D(0, 0,radius)],
               [IMP.algebra.Vector3D(-radius, 0,0), IMP.algebra.Vector3D(0, 0,radius)],
               [IMP.algebra.Vector3D(radius, 0, 0), IMP.algebra.Vector3D(0, 0,-radius)],
               [IMP.algebra.Vector3D(-radius, 0,0), IMP.algebra.Vector3D(0, 0,-radius)])
        if f:
          for d in zip(ds, sites):
            IMP.npctransport.add_sites(f, IMP.core.Typed(d[0]).get_type(),
                                       .5*radius, d[1])
        rs=[]
        for p in [(0,1), (0,2), (0,3), (1,2), (1,3), (2,3)]:
          ps= IMP.npctransport.SitesPairScore(site_range, site_k, nonspec_range,
                                              nonspec_k, soft_sphere_k, sites[p[0]], sites[p[1]])
          ps.set_log_level(IMP.VERBOSE)
          r=IMP.core.PairRestraint(ps, (ds[p[0]], ds[p[1]]))
          rs.append(r)
          r.evaluate(False)
          d=r.create_current_decomposition()
          if d:
            d.set_was_used(False)
        return rs
    def _test_three(self, site_range, site_k, nonspec_range, nonspec_k,
                  soft_sphere_k, dt):
        m= IMP.Model()
        m.set_log_level(IMP.SILENT)
        ds= [self._create_particle(m) for i in range(0,4)]
        ds[0].set_coordinates(IMP.algebra.Vector3D(0,0,0))
        ds[1].set_coordinates(IMP.algebra.Vector3D(2*radius,0,0))
        ds[2].set_coordinates(IMP.algebra.Vector3D(0,0,2*radius))
        ds[3].set_coordinates(IMP.algebra.Vector3D(2*radius,0,2*radius))
        types=[IMP.core.ParticleType(d.get_name()+" type") for d in ds]
        for d in zip(types, ds):
          IMP.core.Typed.setup_particle(d[1], d[0])
        f= RMF.create_rmf_file(self.get_tmp_file_name("glue3.rmf"))
        rs= self._create_restraint_three(m, ds,  site_range, site_k, nonspec_range, nonspec_k,
                                          soft_sphere_k, f)
        bd= IMP.atom.BrownianDynamics(m)
        bd.set_scoring_function(rs)
        bd.set_maximum_time_step(dt)
        w= IMP.npctransport.add_hierarchies(f, ds)
        IMP.rmf.add_restraints(f, rs)
        sos= IMP.rmf.SaveOptimizerState(f)
        bd.add_optimizer_state(sos)
        bd.optimize(3000)
    def test_three(self):
        """Check two interactions"""
        k=500
        dt=IMP.npctransport.get_time_step(1, k, radius)
        self._test_three(radius, k, .2*radius, .5*k, k, dt)
    def _rescore_three(self):
      f= RMF.open_rmf_file_read_only(self.get_tmp_file_name("glue3.rmf"))
      m= IMP.Model()
      ds= IMP.rmf.create_hierarchies(f, m)
      k=500
      rs= self._create_restraint_three(m, ds,  radius, k, .2*radius, .5*k, k)
      IMP.rmf.load_frame(f, 0)
      print rs[0].evaluate(True)
      IMP.rmf.load_frame(f, 1000)
      print rs[0].evaluate(True)
if __name__ == '__main__':
    IMP.test.main()
