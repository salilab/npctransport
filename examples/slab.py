from __future__ import print_function
import IMP.npctransport
import IMP
import RMF
import IMP.rmf
import IMP.atom

height=20
radius=20
k=100

m= IMP.Model()
p= IMP.Particle(m)
IMP.atom.Hierarchy.setup_particle(p)
d=IMP.core.XYZR.setup_particle(p)
d.set_radius(10)
d.set_coordinates_are_optimized(True)
IMP.atom.Diffusion.setup_particle(p)
IMP.atom.Mass.setup_particle(p, 1)
slabss= IMP.npctransport.SlabSingletonScore(height, radius, k)
r= IMP.core.SingletonRestraint(slabss, p, "slab")
bb= IMP.algebra.BoundingBox3D(IMP.algebra.Vector3D(-100, -100, -100),
                              IMP.algebra.Vector3D(100, 100, 100))
bbss= IMP.core.BoundingBox3DSingletonScore(IMP.core.HarmonicUpperBound(0,10),
                                         bb)
bbr= IMP.core.SingletonRestraint(bbss, p, "bb")

nm=IMP.create_temporary_file_name("slab", ".rmf")
rmf= RMF.create_rmf_file(nm)
IMP.rmf.add_hierarchies(rmf, [p])
bbg= IMP.display.BoundingBoxGeometry(bb)
sg= IMP.npctransport.SlabWireGeometry(height, radius, 100)
sg.set_was_used(True);
sgs= sg.get_components()
# silliness we have to do for now
IMP.rmf.add_static_geometries(rmf, [bbg]+sgs)
IMP.rmf.add_restraints(rmf, [bbr, r])

os= IMP.rmf.SaveOptimizerState(m, rmf)

bd=IMP.atom.BrownianDynamics(m)
bd.set_scoring_function([bbr, r])
bd.set_log_level(IMP.SILENT)
bd.add_optimizer_state(os)
bd.set_maximum_time_step(2000)

bd.optimize(100)
print("File is", nm)
