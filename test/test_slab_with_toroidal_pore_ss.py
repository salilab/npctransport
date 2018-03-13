from __future__ import print_function
import IMP
import IMP.test
import IMP.algebra
import IMP.core
import IMP.npctransport
import math
import IMP.display
import random

radius=1
#random.uniform(1,12)
#random.uniform(radius+2, radius+15)
slab_thickness=5
r_h2v_ratio=2 # ration between horizontal minor radius and vertical minor radius (semiaxes) - 1.0 is a ring toroid, >1.0 is a circular torus made of a horizontally elongated ellipse, and <1.0 is made of a vertically elongated ellipse
rv=0.5*slab_thickness # vertical minor radius (semi-axis)
rh=rv*r_h2v_ratio # horizontal minor radius (semi-axis)
slab_radius=5+rh
rv2=rv**2
rh2=rh**2
#random.uniform(5,30)
boxw= 2*max([1.1*slab_radius,slab_thickness])

##
def get_surface_distance_from_axis_aligned_ellipsoid(sphere, origin, rv, rh):
    '''
    return the distance between the surface of specified sphere and
    the surface of an axis aligned ellipsoid with vertical semi-axis rv and
    horizontal semi-axis rh, centered at origin.

    If the sphere penetrates the ellipsoid by some distance d, returns -d.
    '''
    EPS=0.000001
    v=sphere.get_center()-origin
    dXY2= v[0]**2+v[1]**2
    dZ2 = v[2]**2
    #    theta = atan(dXY/(dZ+EPS))
    dv2= dXY2 + dZ2 + EPS
    sinTheta2=dXY2/dv2
    cosTheta2=dZ2/dv2
    cur_r = math.sqrt(rv**2*cosTheta2 + rh**2*sinTheta2)
    dv=math.sqrt(dv2)
    return dv-cur_r-sphere.get_radius()

##
def out_slab(p, ALLOWED_OVERLAP=0.0):
    '''
    verify that spherical (XYZR) particle p is out of slab with
    thickness slab_thickness, and toroidal pore with major radius
    slab_radius and minor radius 0.5*slab_thickness.

    penetration of ALLOWED_OVERLAP is allowed (would still result in True). Use negative value to force extra slack
    '''
    EPS=0.000001
    R=slab_radius # major radius of torus
    d= IMP.core.XYZR(p)
    print("out_slab begins: ", d,R,rv, rh)
    xyz= d.get_coordinates()
    if xyz[2] - d.get_radius() + ALLOWED_OVERLAP > rv:
        print("Above slab")
        return True
    if xyz[2] + d.get_radius() - ALLOWED_OVERLAP < -rv:
        print("Below slab")
        return True
    print("In slab vertically (but still might be in pore)")
    xy0= IMP.algebra.Vector3D(xyz[0], xyz[1], 0.0)
    rxy0= xy0.get_magnitude()
    print("xy0", xy0, "magnitude", rxy0)
    if rxy0 + d.get_radius() > R:
        print("Sphere overlaps slab outside the torus major radius")
        return False
    if(rxy0>EPS):
        xy0_major= xy0 * ( R / rxy0 ) # same direction as xy0 but on major radius
    else:
        xy0_major=IMP.algebra.Vector3D(0, R, 0.0)
    print("xy0_major", xy0_major)
    distance = get_surface_distance_from_axis_aligned_ellipsoid \
               (d.get_sphere(), xy0_major, rv, rh)
    print("Sphere surface distance to torus", distance)
    is_overlap = (distance + ALLOWED_OVERLAP < 0)
    print("Is overlapping slab somewhere within torus major radius? ", is_overlap)
    return not is_overlap

##
class ConeTests(IMP.test.TestCase):
    def test_slab_singleton_score(self):
        """Check slab singleton score"""
        IMP.set_log_level(IMP.SILENT)
        print("radius", radius, "slab radius", slab_radius, "slab_thickness", slab_thickness)
        m= IMP.Model()
        m.set_log_level(IMP.SILENT)
        p= IMP.Particle(m)
        d= IMP.core.XYZR.setup_particle(p)
        d.set_coordinates_are_optimized(True)
        d.set_radius(radius)
        bb= IMP.algebra.BoundingBox3D(0.5*IMP.algebra.Vector3D(-boxw, -boxw, -boxw),
                                      0.5*IMP.algebra.Vector3D(boxw,boxw,boxw))
        p_slab= IMP.Particle(m, "slab")
        IMP.npctransport.SlabWithToroidalPore.setup_particle \
              (p_slab, slab_thickness, slab_radius, r_h2v_ratio)
        self.assert_(IMP.npctransport.SlabWithToroidalPore.get_is_setup(p_slab))
        # test cast to slab with pore
        slab= IMP.npctransport.SlabWithToroidalPore(p_slab)
        self.assertEqual(slab.get_pore_radius(),
                         slab_radius)
        self.assertEqual(slab.get_thickness(),
                         slab_thickness);
        self.assertEqual(slab.get_minor_radius_h2v_aspect_ratio(),
                         r_h2v_ratio)
        self.assertEqual(slab.get_vertical_minor_radius(),
                         rv);
        self.assertEqual(slab.get_horizontal_minor_radius(),
                         rh);
        slabps= IMP.npctransport.SlabWithToroidalPorePairScore \
                (1.0)
        slabps.set_log_level(IMP.SILENT)
        r= IMP.core.PairRestraint(m,
                                  slabps,
                                  [p_slab.get_index(),p.get_index()],
                                  "slab")
        while out_slab(p, ALLOWED_OVERLAP=0.5):
            d.set_coordinates(IMP.algebra.get_random_vector_in(bb))
        print("Penetrating slab: ", d.get_coordinates())
        pym_fname="tmp.pym" #self.get_tmp_file_name("slabps.pym")
        w= IMP.display.PymolWriter(pym_fname)
        w.set_frame(0)
        g=IMP.core.XYZRGeometry(d)
        sg= IMP.npctransport.SlabWithToroidalPoreWireGeometry(slab_thickness,
                                                              slab_radius,
                                                              rh,
                                                              boxw)
        w.add_geometry([g, sg])
        cg= IMP.core.SteepestDescent(m)
        cg.set_scoring_function(r)
        cg.set_log_level(IMP.SILENT)
        f=1
        for i in range(0,100):
            s=cg.optimize(10)
            w.set_frame(f)
            f=f+1
            w.add_geometry([g, sg])
            print(i, d.get_coordinates())
            if s==0:
                self.assert_(out_slab(d, ALLOWED_OVERLAP=0.0))
                break
            else:
                self.assert_(not out_slab(d, ALLOWED_OVERLAP=0.0))
        print(d.get_coordinates())

if __name__ == '__main__':
    IMP.test.main()
