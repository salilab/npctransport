#!/usr/bin/env python
"""Add cylinders of a given radius to the passed file."""
import IMP
import RMF
import IMP.rmf
from IMP.npctransport import *
import sys
import IMP.display
import numpy as np

def clone_rmf_static(in_name,out_name):
    return in_file, out_file


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
    '''
    node - rmf node to scan
    cf - cylinder factory
    cdf - colored factory
    tf - typed factory
    types - list of type names for which to apply method
    '''
    children = node.get_children()
    ret = []
    #print "inspecting", node.get_name()
    if len(children)==0:
        return ret
    if has_depth_with_site(node, 3) and tf.get_is(children[0]):
        tf_type = tf.get(children[0]).get_type_name()
        if (tf_type in types):
            #            print tf_type, "is of right type"
            for i in range(0, len(children) - 1):
                cyl = node.add_child("cylinder", RMF.GEOMETRY)
                cdf.get(cyl).set_static_rgb_color(RMF.Vector3(color[0],color[1],color[2]))
                ret.append((cyl, children[i], children[i+1]))
                #                print "adding for", ret[-1], "depth", depth, "d3under?", has_depth_with_site(node, 3)
                cf.get(cyl).set_radius(radius)
    for c in children:
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
    coords_list= [RMF.Vector3(cep0), RMF.Vector3(cep1)]
#    coords=[[cep0[0], cep1[0]], [cep0[1], cep1[1]], [cep0[2], cep1[2]]]
    cf.get(nh).set_frame_coordinates_list(coords_list)


def _smooth(node, tf, rff, xyz_dict, n=10, is_write=True):
    '''
    node - node being scanned
    tf - typed factory
    rff - reference frame factory
#    bf - ball factory
    xyz_dict - dictinory from node ids to list of up to n xyz coordinates
    n - smoothing window size
    '''
    if tf.get_is(node) and rff.get_is(node):
        node_id= node.get_id()
        rf= rff.get(node)
        xyz= np.array([t for t in rf.get_frame_translation()])
        if not node_id in xyz_dict:
            xyz_dict[node_id]= [xyz]
        else:
            xyz_dict[node_id].append(xyz)
        if len(xyz_dict[node_id])>n:
            xyz_dict[node_id]= xyz_dict[node_id][-n:]
        smooth_xyz= np.mean(xyz_dict[node_id],axis=0)
        if is_write:
            # print "xyz", xyz
            # print "XYZ_LIST", xyz_dict[node_id]
            # print "smooth_xyz", smooth_xyz
            # print len(smooth_xyz)
            rf.set_frame_translation(
                RMF.Vector3(smooth_xyz[0],
                            smooth_xyz[1],
                            smooth_xyz[2]))
    for c in node.get_children():
        _smooth(c, tf, rff, xyz_dict, n)


def _recolor(node, tf, cf, types, color):
    children= node.get_children()
    if tf.get_is(node):
        if tf.get(node).get_type_name() in types:
            cd= cf.get(node)
            cd.set_static_rgb_color(RMF.Vector3(color[0],color[1],color[2]))
    for c in children:
#        print "recolor",c,color
        _recolor(c, tf, cf, types, color)

def _resize_sites(node, bf, nr):
    children = node.get_children()
    if node.get_name() == "site" and bf.get_is(node):
        d= bf.get(node)
        d.set_static_radius(nr)
    for c in children:
        _resize_sites(c, bf, nr)

def _get_fg_and_floater_types(ref_output):
    fg_types = []
    kap_types= []
    inert_types= []
    output = None
    try:
        if(ref_output <> ""):
            output=Output()
            FILE=open(ref_output,"rb")
            output.ParseFromString(FILE.read())
    except:
        print "Couldn't read '" + ref_output + "'"
        raise
    if(output == None):
        for i in range(0, get_number_of_types_of_fg()):
            fg_types.append( get_type_of_fg(i).get_string() )
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
        #        for i in range(0, get_number_of_types_of_float()):
        #           floater_types.append( get_type_of_float(i).get_string() )
    else:
        a = output.assignment
        for fg in a.fgs:
            fg_types.append(fg.type)
            print "Added fg type", fg.type
        for floater in a.floaters:
            if(floater.interactions.value>0):
                kap_types.append(floater.type)
            else:
                inert_types.append(floater.type)
    print "FGs:", fg_types
    print "Kaps:", kap_types
    print "Interts:", inert_types
    return fg_types, kap_types, inert_types


fg_color=[255.0/255, 204.0/255 , 102.0/255]
kap_color=[220.0/255, 0, 64.0/255]
inert_color=[70.0/255, 80.0/255, 220.0/255]
def main():
    IMP.add_string_flag("input_rmf", "", "The input RMF file.")
    IMP.add_string_flag("output_rmf", "", "The output RMF file in which to add cylinders.")
    IMP.add_string_flag("ref_output", "", "reference output file from which info e.g. fg nup types can be extracted")
    IMP.add_float_flag("radius", 5, "The radius of the cylinder.")
    IMP.add_float_flag("site_radius", 2, "The radius of the sites.")
    IMP.add_bool_flag("recolor_fgs", "recolor fg nup chains")
    IMP.add_bool_flag("recolor_floats", "recolor floating (diffusing) molecules")
    IMP.add_int_flag("smooth_n_frames", 1,
                     "smooth by averaging over a window of specified frames (but retaining the same number of expected frames")
    IMP.add_int_flag("skip_n_frames", 1,
                     "skip every n frames of output")
    IMP.setup_from_argv(sys.argv, "Prettify a movie")
    in_fname= IMP.get_string_flag("input_rmf")
    out_fname= IMP.get_string_flag("output_rmf")
    ref_output = IMP.get_string_flag("ref_output")
    radius= IMP.get_float_flag("radius")
    # Prepare out file with same static info as in file:
    in_fh = RMF.open_rmf_file_read_only(in_fname)
    out_fh = RMF.create_rmf_file(out_fname)
    print("creating file", out_fname)
    RMF.clone_file_info(in_fh, out_fh)
    RMF.clone_hierarchy(in_fh, out_fh)
    RMF.clone_static_frame(in_fh, out_fh)
    print "opened", in_fh.get_name()
    cf = RMF.CylinderFactory(out_fh)
    rff = RMF.ReferenceFrameFactory(out_fh)
    tf = RMF.TypedFactory(out_fh)
    bf = RMF.BallFactory(out_fh)
    cdf = RMF.ColoredFactory(out_fh)
    ipf = RMF.IntermediateParticleFactory(out_fh)
    fg_types, kap_types, inert_types = _get_fg_and_floater_types( ref_output )
#    out_fh.set_current_frame(RMF.ALL_FRAMES)
    # Modify static information:
    if(IMP.get_bool_flag("recolor_floats")):
        _recolor(out_fh.get_root_node(), tf, cdf, kap_types, kap_color)
        _recolor(out_fh.get_root_node(), tf, cdf, inert_types, inert_color)

    _resize_sites(out_fh.get_root_node(), bf, IMP.get_float_flag("site_radius"))
    cylinders = []
    for i, type in enumerate(fg_types):
        color = IMP.display.get_display_color(i)
        rgb = [color.get_red(), color.get_green(), color.get_blue()]
        cylinders += _add_nodes(out_fh.get_root_node(), cf, cdf, tf, [type], radius, rgb) # fg_color
        if(IMP.get_bool_flag("recolor_fgs")):
            #            print "Recoloring",type, rgb
            _recolor(out_fh.get_root_node(), tf, cdf, [type], rgb) # fg_color
    # Clone and modify per-frame information:
    smooth_xyz_dict={}
    smooth_n_frames=IMP.get_int_flag("smooth_n_frames")
    skip_n_frames=IMP.get_int_flag("skip_n_frames")
    print "Skip interval:", skip_n_frames, "frames"
    for f_id, f in enumerate(in_fh.get_frames()):
        is_write= f_id % skip_n_frames == 0
        in_fh.set_current_frame(f)
        if not is_write:
            if(smooth_n_frames>1):
                _smooth(in_fh.get_root_node(),
                        tf,
                        rff,
                        smooth_xyz_dict,
                        n=smooth_n_frames,
                        is_write=False)
            print("skipping frame", f, f_id)
            continue
        print("cloning frame", f, f_id)
        out_fh.add_frame(in_fh.get_name(f), in_fh.get_type(f))
        RMF.clone_loaded_frame(in_fh, out_fh)
        if(smooth_n_frames>1):
            _smooth(out_fh.get_root_node(),
                    tf,
                    rff,
                    smooth_xyz_dict,
                    n=smooth_n_frames,
                    is_write=True)
        for c in cylinders:
            _set_cylinder(c, cf, rff)
        DEBUG=False
        if DEBUG and f_id >= 5*skip_n_frames:
            break

main()
