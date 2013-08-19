#!/usr/bin/env python
"""Add cylinders of a given radius to the passed file."""
import IMP.base
import RMF
import IMP.npctransport
import sys

def has_depth_with_site(root, i):
  """ returns true if node subtree thru first child is at least i
      levels, including the root node itself, and the lead is a site """
#  print root, i, len(root.get_children())
  if (i==1) and root.get_name()=="site":
    return True
  c = root.get_children()
  if len(c) == 0:
    return False
  return has_depth_with_site(c[0], i-1)


def _add_nodes(node, cf, cdf, tf, types, radius, color, depth=0):
  children = node.get_children()
  ret = []
  #print "inspecting", node.get_name()
  if len(children)==0:
    return ret
  #if has_depth_with_site(node, 3) and
  if tf.get_is(children[0]):
    if (tf.get(children[0]).get_type_name() in types):
      print "it is"
      for i in range(0, len(children) - 1):
        cyl = node.add_child("cylinder", RMF.GEOMETRY)
        cdf.get(cyl).set_rgb_color(color)
        ret.append((cyl, children[i], children[i+1]))
        print "adding for", ret[-1], "depth", depth, "d3under?", has_depth_with_site(node, 3)
        cf.get(cyl).set_radius(radius)
        return ret # only one cylinder
#      if(depth < 2):
  for c in children:
    if(len(ret)>0):
      return ret
    ret += _add_nodes(c, cf, cdf, tf, types, radius, color, depth+1)
  return ret

def _set_cylinder(cylinder_descriptor, cf, ipf):
  """
  draws a cylinder between two nodes with refframes
  cylinder descriptor - tuple (cylinder_node, node1, node2)
  cf - cylinder factory
  ipf - refframe factory
  """
  nh = cylinder_descriptor[0]
  ep0 = cylinder_descriptor[1]
  ep1 = cylinder_descriptor[2]
  cep0 = ipf.get(ep0).get_translation()
  cep1 = ipf.get(ep1).get_translation()
  coords=[[cep0[0], cep1[0]], [cep0[1], cep1[1]], [cep0[2], cep1[2]]]
  cf.get(nh).set_coordinates(coords)


def _recolor(node, tf, cf, types, color):
  children= node.get_children()
  if tf.get_is(node):
    if tf.get(node).get_type_name() in types:
      cd= cf.get(node)
      cd.set_rgb_color(color)
  for c in children:
    _recolor(c, tf, cf, types, color)

def _resize_sites(node, bf, nr):
  children = node.get_children()
  if node.get_name() == "site" and bf.get_is(node):
    d= bf.get(node)
    d.set_radius(nr)
  for c in children:
    _resize_sites(c, bf, nr)

fg_color=[255.0/255, 204.0/255 , 102.0/255]
kap_color=[128.0/255, 0, 64.0/255]
crap_color=[70.0/255, 90.0/255, 220.0/255]
def main():
  IMP.base.add_string_flag("input", "", "The RMF file to add cylinders to.")
  IMP.base.add_float_flag("radius", 5, "The radius of the cylinder.")
  IMP.base.add_float_flag("site_radius", 2, "The radius of the sites.")
  IMP.base.add_bool_flag("recolor_fgs", "recolor fg nup chains")
  IMP.base.add_bool_flag("recolor_floats", "recolor floating (diffusing) molecules")
  IMP.base.setup_from_argv(sys.argv, "Prettify a movie")
  fh= RMF.open_rmf_file(IMP.base.get_string_flag("input"))
  radius= IMP.base.get_float_flag("radius")
  print "opened", fh.get_name()
  cf = RMF.CylinderFactory(fh)
  rff = RMF.ReferenceFrameFactory(fh)
  tf = RMF.TypedFactory(fh)
  bf = RMF.BallFactory(fh)
  cdf = RMF.ColoredFactory(fh)
  ipf = RMF.IntermediateParticleFactory(fh)
  fg_types = [IMP.npctransport.get_type_of_fg(i).get_string() for i in range(0, IMP.npctransport.get_number_of_types_of_fg())]
  fg_types = fg_types + [ "Nup57_16copies_chimera",
                          "Nup49_16copies",
                          "Nsp1_16copies_1",
                          "Nsp1_16copies_2",
                          "Nup159_8copies",
                          "Nup116_8copies_chimera",
                          "Nup42_8copies_chimera",
                          "Nup100_8copies_chimera",
                          "Nup145N_8copies_1_chimera",
                          "Nup145N_8copies_2_chimera",
                          "Nup1_8copies",
                          "Nup60_8copies"]
  print fg_types
  float_types = [IMP.npctransport.get_type_of_float(i).get_string() for i in range(0, IMP.npctransport.get_number_of_types_of_float())]
  fh.set_current_frame(RMF.ALL_FRAMES)

#  if(IMP.base.get_bool_flag("recolor_fgs")):
#    _recolor(fh.get_root_node(), tf, cdf, fg_types, fg_color)
#  if(IMP.base.get_bool_flag("recolor_floats")):
#    _recolor(fh.get_root_node(), tf, cdf, float_types[0:1], kap_color)
#    _recolor(fh.get_root_node(), tf, cdf, float_types[1:], crap_color)

#  _resize_sites(fh.get_root_node(), bf, IMP.base.get_float_flag("site_radius"))

  cylinders = _add_nodes(fh.get_root_node(), cf, cdf, tf, fg_types, radius, fg_color)
  for i in range(0, fh.get_number_of_frames()):
    print "frame", i
    fh.set_current_frame( RMF.FrameID(i) )
    for c in cylinders:
      _set_cylinder(c, cf, rff)

main()