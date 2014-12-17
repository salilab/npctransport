from __future__ import print_function
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
    def _randomize(self, rbs, sites, bb):
        '''
        randomize reference frames of rbs such that each rb and
        its first site are not too far or too close to all other
        rbs or their first sites.
        '''
        for i in range(len(rbs)):
            ok=False
            failures=0
            while not ok:
                tr= IMP.algebra.get_random_vector_in(bb)
                r= IMP.algebra.get_random_rotation_3d()
                trans= IMP.algebra.Transformation3D(r, tr)
                rbs[i].set_reference_frame(IMP.algebra.ReferenceFrame3D(trans))
                ok=True
                for orb, s in zip(rbs[0:i], sites[0:i]):
                    if failures>200:
                        return self._randomize(rbs, sites, bb) # retry
                    d= IMP.core.get_distance(IMP.core.XYZR(rbs[i]),
                                             IMP.core.XYZR(orb))
                    if d <0 or d >radius:
                        failures+=1
                        ok=False
                        break
                    sp=rbs[i].get_reference_frame().get_global_coordinates(sites[i][0])
                    spo= orb.get_reference_frame().get_global_coordinates(s[0])
                    ds= IMP.algebra.get_distance(sp, spo)
                    print("i=",i, "d=", d, "dsites=", ds)
                    if ds > (radius/2.0):
                        failures+=1
                        ok=False

    def _show(self, rbs, sites, w):
        for i,r in enumerate(rbs):
            g= IMP.npctransport.SitesGeometry(r, sites[i])
            w.add_geometry(g)

    def test_sites_pair_score(self):
        """Check sites pair score"""
        print("Check sites pair score")
        m= IMP.Model()
        m.set_log_level(IMP.SILENT)
        rb0= self._create_particle(m)
        rb1= self._create_particle(m)
        bb= IMP.algebra.get_cube_3d(15)
        s0=[IMP.algebra.Vector3D(radius, 0,0)]
        s1=[IMP.algebra.Vector3D(0, radius,0)]
        r_sites = 1000
        k_sites = 10
        r_nonspec_atr = 0
        k_nonspec_atr = 0
        k_rep = 10
        ps= IMP.npctransport.SitesPairScore(r_sites,k_sites,
                                            r_nonspec_atr, k_nonspec_atr, k_rep,
                                            s0,
                                            s1)
        IMP.set_log_level(IMP.SILENT)
        r= IMP.core.PairRestraint(ps, (rb0, rb1))
        m.add_restraint(r)
        w= IMP.display.PymolWriter("out.pym")
        w.set_frame(0)
        self._randomize([rb0, rb1], [s0, s1], bb)
        self._show([rb0, rb1], [s0, s1], w)
        cg= IMP.core.ConjugateGradients(m)
        print("SCORE=",m.evaluate(True))

        rb0.get_particle().show()
        rb1.get_particle().show()

        cg.optimize(1000)
        print("SCORE=",m.evaluate(True))
        rb0.get_particle().show()
        rb1.get_particle().show()
        w.set_frame(1)
        self._show([rb0, rb1], [s0, s1], w)
        rf0= rb0.get_reference_frame()
        rf1= rb1.get_reference_frame()
        s0g= rf0.get_global_coordinates(s0[0])
        s1g= rf1.get_global_coordinates(s1[0])
        d= IMP.algebra.get_distance(s0g, s1g)
        print("Sites", s0g, s1g, "d=", d)
        self.assert_(d < 1)
if __name__ == '__main__':
    IMP.test.main()
