import RMF
import sys
import os.path

for fn in sys.argv[1:]:
    f = RMF.open_rmf_file_read_only(fn)
    fu = RMF.create_rmf_file(os.path.splitext(fn)[0] + ".updated.rmf")
    RMF.clone_file_info(f, fu)
    fu.set_producer("missing particle updated\n")
    RMF.clone_hierarchy(f, fu)
    RMF.clone_static_frame(f, fu)

    pcat = fu.get_category("physics")
    rkey = fu.get_key(pcat, "radius", RMF.float_traits)
    ckey = fu.get_key(pcat, "coordinates", RMF.vector3_traits)
    for nid in fu.get_node_ids():
        n = fu.get_node(nid)
        if n.get_value(rkey) is not None:
            n.set_static_value(ckey, RMF.Vector3(0, 0, 0))

    for fr in f.get_frames():
        f.set_current_frame(fr)
        fu.add_frame(f.get_current_frame_name(), f.get_current_frame_type())
        RMF.clone_loaded_frame(f, fu)
