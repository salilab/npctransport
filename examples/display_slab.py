import IMP.npctransport
import IMP
import RMF
import IMP.rmf
import IMP.atom
import IMP.display

height=20
radius=30
particle_radius=10
sz=60
use_rmf=True

m= IMP.Model()
m.set_log_level(IMP.WARNING)
p= IMP.Particle(m)
d=IMP.core.XYZR.setup_particle(p)
d.set_radius(10)
d.set_coordinates_are_optimized(True)
IMP.atom.Mass.setup_particle(p, 1)
IMP.atom.Hierarchy.setup_particle(p)

slabss= IMP.npctransport.SlabSingletonScore(height, radius, 1)
r= IMP.core.SingletonRestraint(slabss, p, "slab")
r.set_log_level(IMP.WARNING)
nm=IMP.base.create_temporary_file_name("display_slab", ".pym")
w= IMP.display.create_writer(nm)
if use_rmf:
    rnm=IMP.base.create_temporary_file_name("display_slab", ".rmf")
    rmf= RMF.create_rmf_file(rnm)
    #IMP.rmf.add_hierarchies(rmf, [p])
sg= IMP.npctransport.SlabWireGeometry(height, radius, 100)
sg.set_was_used(True);
if use_rmf:
    sg.set_color(IMP.display.Color(0,1,0))
    sgs= sg.get_components()
    nh= rmf.get_root_node().add_child("slab", RMF.GEOMETRY)
    IMP.rmf.add_geometries(nh, sgs)
    #[sc.show() for sc in sgs]
w.add_geometry(sg)

bb=IMP.algebra.BoundingBox3D(IMP.algebra.Vector3D(0, 0, 0),
                             IMP.algebra.Vector3D(sz, sz, sz))

#nh= rmf.get_root_node().add_child("derivatives", RMF.GEOMETRY)

gs=[]
if use_rmf:
    nh= rmf.get_root_node().add_child("deriv", RMF.GEOMETRY)
for c in IMP.algebra.get_grid_interior_cover_by_spacing(bb, 5):
    #print c
    d.set_coordinates(c)
    v=r.evaluate(True)
    if v ==0:
        continue
    color= IMP.display.get_hot_color(v/20.0)
    seg= IMP.algebra.Segment3D(c, c+d.get_derivatives())
    #print c, seg
    g= IMP.display.SegmentGeometry(seg)
    g.set_color(color)
    g.set_name("deriv")
    gs.append(g)
if use_rmf:
    IMP.rmf.add_geometries(nh, gs)
w.add_geometry(gs)
#bb=IMP.algebra.BoundingBox3D(IMP.algebra.Vector3D(-sz, -sz, -sz),
#                             IMP.algebra.Vector3D(0, 0, 0))
if use_rmf:
    nh= rmf.get_root_node().add_child("displacements", RMF.GEOMETRY)

gs=[]
for c in IMP.algebra.get_grid_interior_cover_by_spacing(bb, 5):
    dv= slabss.get_displacement_direction(c)
    #.get_unit_vector()
    dm= slabss.get_displacement_magnitude(c)
    #print c, dm, dv
    if dm > particle_radius:
        continue
    color= IMP.display.get_hot_color(-(dm-particle_radius)/(height+particle_radius))
    m.evaluate(True)
    seg= IMP.algebra.Segment3D(c, c+dv)
    g= IMP.display.SegmentGeometry(seg)
    g.set_color(color)
    g.set_name("displace")
    gs.append(g)
if use_rmf:
    IMP.rmf.add_geometries(nh, gs)
w.add_geometry(gs)

if use_rmf:
    IMP.rmf.save_frame(rmf, 0)
    print "chimera", rnm
print "pymol", nm
