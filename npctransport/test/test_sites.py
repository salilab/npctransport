import IMP
import IMP.test
import IMP.npctransport
import math

radius=5

class ConeTests(IMP.test.TestCase):
    def _create_particle(self, m):
        p= IMP.Particle(m)
        d= IMP.core.XYZR.setup_particle(p)
        d.set_radius(radius)
        rb= IMP.core.RigidBody.setup_particle(p, IMP.algebra.ReferenceFrame3D())
        rb.set_coordinates_are_optimized(True)
        return rb
    def _randomize(self, rps, sites, bb):
        for i in range(len(rps)):
            ok=False
            failures=0
            while not ok:
                tr= IMP.algebra.get_random_vector_in(bb)
                r= IMP.algebra.get_random_rotation_3d()
                trans= IMP.algebra.Transformation3D(r, tr)
                rps[i].set_reference_frame(IMP.algebra.ReferenceFrame3D(trans))
                ok=True
                for orb, s in zip(rps[0:i], sites[0:i]):
                    if failures>200:
                        raise RuntimeError("too many failures with "+str(IMP.core.XYZR(rps[i])) + " and "+str(IMP.core.XYZR(orb)))
                    d= IMP.core.get_distance(IMP.core.XYZR(rps[i]),
                                             IMP.core.XYZR(orb))
                    if d <0 or d >radius:
                        ok=False
                        failures+=1
                        break
                    sp=rps[i].get_reference_frame().get_global_coordinates(sites[i][0])
                    spo= orb.get_reference_frame().get_global_coordinates(s[0])
                    ds= IMP.algebra.get_distance(sp, spo)
                    print d, ds
                    if ds >5:
                        failures+=1
                        ok=False
    def _show(self, rps, sites, w):
        for i,r in enumerate(rps):
            g= IMP.npctransport.SitesGeometry(r, sites[i])
            w.add_geometry(g)
    def test_cone_construction(self):
        """Check sites pair score"""
        m= IMP.Model()
        m.set_log_level(IMP.SILENT)
        rp0= self._create_particle(m)
        rp1= self._create_particle(m)
        bb= IMP.algebra.get_cube_3d(15)
        s0=[IMP.algebra.Vector3D(radius, 0,0)]
        s1=[IMP.algebra.Vector3D(0, radius,0)]
        ps= IMP.npctransport.SitesPairScore(1000, 10,
                                            0,0,10,
                                            s0,
                                            s1)
        ps.set_log_level(IMP.VERBOSE)
        r= IMP.core.PairRestraint(ps, (rp0, rp1))
        m.add_restraint(r)
        w= IMP.display.PymolWriter("out.pym")
        w.set_frame(0)
        self._randomize([rp0, rp1], [s0, s1], bb)
        self._show([rp0, rp1], [s0, s1], w)
        cg= IMP.core.ConjugateGradients(m)

        m.evaluate(True)
        rp0.get_particle().show()
        rp1.get_particle().show()

        cg.optimize(1000000)
        rp0.get_particle().show()
        rp1.get_particle().show()
        w.set_frame(1)
        self._show([rp0, rp1], [s0, s1], w)
        rf0= rp0.get_reference_frame()
        rf1= rp1.get_reference_frame()
        s0g= rf0.get_global_coordinates(s0[0])
        s1g= rf1.get_global_coordinates(s1[0])
        d= IMP.algebra.get_distance(s0g, s1g)
        print s0g, s1g, d
        self.assert_(d < 1)
    def _test_cone_construction(self):
        """Check sites pair score stable"""
        m= IMP.Model()
        rp0= self._create_particle(m)
        rp1= self._create_particle(m)
        bb= IMP.algebra.get_cube_3d(10)
        rp0.set_coordinates(IMP.algebra.Vector3D(-radius, 0,0))
        rp1.set_coordinates(IMP.algebra.Vector3D(radius, 0,0))
        s0=[IMP.algebra.Vector3D(radius, 0,0)]
        s1=[IMP.algebra.Vector3D(-radius, 0,0)]
        ps= IMP.npctransport.SitesPairScore(1000, 10,
                                             s0,
                                             s1)
        r= IMP.core.PairRestraint(ps, (rp0, rp1))
        m.add_restraint(r)
        w= IMP.display.PymolWriter("outs.pym")
        w.set_frame(0)
        self._show([rp0, rp1], [s0, s1], w)
        cg= IMP.core.ConjugateGradients(m)
        cg.optimize(1000)
        w.set_frame(1)
        self._show([rp0, rp1], [s0, s1], w)
        rf0= rp0.get_reference_frame()
        rf1= rp1.get_reference_frame()
        s0g= rf0.get_global_coordinates(s0[0])
        s1g= rf1.get_global_coordinates(s1[0])
        d= IMP.algebra.get_distance(s0g, s1g)
        print s0g, s1g, d
        self.assert_(d < 1)
if __name__ == '__main__':
    IMP.test.main()
