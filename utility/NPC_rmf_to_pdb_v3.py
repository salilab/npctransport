#!/usr/bin/env python
"""Add cylinders of a given radius to the passed file."""
from IMP.npctransport import *
import IMP
import RMF
import sys
import IMP.display
import numpy as np
import os

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


def _add_nodes(node, tf, type_prefixes, depth=0):
    '''
    node - rmf node to scan
    tf - typed factory
    type_prefixes - list of full type prefixes (e.g. "Nup1" for "Nup1N")

    adds only nodes whose type name begins with any of the specified type prefixes
    '''
    children = node.get_children()
    ret = []
    #print "inspecting", node.get_name()
    if len(children)==0:
        return ret
    if has_depth_with_site(node, 3) and tf.get_is(children[0]):
        child_type = tf.get(children[0]).get_type_name()
        if any([child_type.startswith(tp) for tp in type_prefixes]):
            ret.append(children)
    for c in children:
        ret += _add_nodes(c, tf,  type_prefixes, depth+1)
    return ret

def _print_atom(node, rf, atom_id, res_id, chain_id, fg_type):
    """
    prints atom entry for specified node as a CA atom
    node - node name
    rf - refframe factory (intermediate particle factory)
    """


def _get_fg_and_floater_types(ref_output):
    fg_types = []
    kap_types= []
    inert_types= []
    output = None
    try:
        if(ref_output != ""):
            output=Output()
            FILE=open(ref_output,"rb")
            output.ParseFromString(FILE.read())
    except:
        print(f"Couldn't read '{ref_output}'")
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
#            print("Added fg type", fg.type)
        for floater in a.floaters:
            if(floater.interactions.value>0):
                kap_types.append(floater.type)
            else:
                inert_types.append(floater.type)
#    print("FGs:", fg_types)
#    print("Kaps:", kap_types)
#    print("Inerts:", inert_types)
    return fg_types, kap_types, inert_types

def my_main():
    fg_color=[255.0/255, 204.0/255 , 102.0/255]
    kap_color=[220.0/255, 0, 64.0/255]
    inert_color=[70.0/255, 80.0/255, 220.0/255]
    IMP.add_string_flag("input_rmf", "", "The input RMF file.")
    IMP.add_string_flag("output_pdb", "", "The output PDB file")
    IMP.add_string_flag("ref_output", "", "reference output file from which info e.g. fg nup types can be extracted")
    IMP.add_int_flag("skip_n_frames", 1,
                     "skip every n frames of output; if 0, include only first frame")
    IMP.add_int_flag("first_frame", 0,
                     "start from speicifed frame")
    IMP.setup_from_argv(sys.argv, "Convert an rmf to a pdb file.")
    in_fname= IMP.get_string_flag("input_rmf")
    out_fname= IMP.get_string_flag("output_pdb")
    ref_output = IMP.get_string_flag("ref_output")
    # Prepare out file with same static info as in file:
    print("Reading", in_fname)
    in_fh = RMF.open_rmf_file_read_only(in_fname)
    rff = RMF.ReferenceFrameFactory(in_fh)
    tf = RMF.TypedFactory(in_fh)
    bf = RMF.BallFactory(in_fh)
    ipf = RMF.IntermediateParticleFactory(in_fh)
    fg_types, kap_types, inert_types = _get_fg_and_floater_types( ref_output )
    type2chains={}
#    print("fg_types", fg_types)
    for i, fg_type in enumerate(fg_types):
        type2chains[fg_type] = _add_nodes(in_fh.get_root_node(), tf, [fg_type])
#        print(type2chains[fg_type])
    skip_n_frames=IMP.get_int_flag("skip_n_frames")
    first_frame=IMP.get_int_flag("first_frame")
    if(skip_n_frames>0):
        print("Skip interval:", skip_n_frames, "frames")
    else:
        print("Showing only first frame")
    for f_id, f in enumerate(in_fh.get_frames()):
        is_write= (f_id>=first_frame) and ((skip_n_frames==0) or (f_id%skip_n_frames == 0))
        in_fh.set_current_frame(f)
        if not is_write:
            print("skipping frame", f, f_id)
            continue
        fname= "%s_f%05d.pdb" % (os.path.splitext(out_fname)[0],f_id)
        F= open(fname, 'w')
        print("COMMENT writing frame", f, f_id)
        chain_id="A"
        atom_id=1
        model_id=1
        spoke_id=1                    
        for fg_type, chains in type2chains.items():
            # print model id to file F
            print(f"MODEL {model_id:d}", file=F)
            for chain in chains:
                if chain_id=='I': # no more than 8 chains per model
                    chain_id='A'
                    print("ENDMDL", file=F)
                    model_id=model_id+1
                    print(f"MODEL {model_id:d}", file=F)
                for i,node in enumerate(chain):
                    res_id=i+1
                    coord = rff.get(node).get_translation()
                    print("ATOM  {:5d}  CA  GLY {:1s}{:4d}    {:8.2f}{:8.2f}{:8.2f}   1.0                         # {:s}".format \
                        (atom_id,  chain_id, res_id,
                           coord[0], coord[1], coord[2],
                           fg_type + f" spoke{spoke_id}"), 
                        file=F)
                    atom_id=atom_id+1
                print("TER", file=F)
                chain_id=chr(ord(chain_id)+1)
                spoke_id = 1 + (spoke_id % 8)
            chain_id='A'
            print("ENDMDL", file=F)
            model_id=model_id+1
        if skip_n_frames==0:
            break
        DEBUG=False
        if (DEBUG and f_id >= 5*skip_n_frames):
            break
        print("Finished frame {:d}".format(f_id))

if __name__ == '__main__':
    my_main()
