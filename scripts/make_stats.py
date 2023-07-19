# #!/netapp/sali/barak/Runs/NPC_BD/Cyl2/tools/imppy.sh /usr/bin/python
from IMP.npctransport import *
import numpy
import math
import sys
import re
import numbers
import os
import glob
import traceback
KEY_CAPTIONS = None

def query_yes_no(question, default=None):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is one of "yes" or "no".
    """
    valid = {"yes":True,   "y":True,  "ye":True,
             "no":False,     "n":False}
    if default == None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "\
                             "(or 'y' or 'n').\n")

def is_crap(typename):
    return re.search("crap", typename) != None

def is_fg(typename):
    return not is_kap(typename) and not is_crap(typename)

def is_kap(typename):
    return re.search("kap", typename) != None

def is_kap_fg_interaction(i):
    """
    i - protobuf interaction

    return true if interaction between kap and fg
    """
    return ( is_fg(i.type0) and is_kap(i.type1) ) or \
        ( is_kap(i.type0) and is_fg(i.type1) )

def is_crap_fg_interaction(i):
    """
    i - protobuf interaction

    return true if interaction between kap and fg
    """
    return ( is_fg(i.type0) and is_crap(i.type1) ) or \
        ( is_crap(i.type0) and is_fg(i.type1) )

def get_interaction_k(ints):
    """
    ints - protobuf Assignment.Interactions

    returns the kap-fg and fg-fg interaction k from interactions
    """
    kap_k = None
    fg_k = None
    for i in ints:
        if is_kap_fg_interaction(i) and kap_k is None:
            kap_k= i.interaction_k.value
        if is_fg(i.type0) and is_fg(i.type1) and fg_k is None:
            fg_k= i.interaction_k.value
        if(kap_k <> None and fg_k <> None):
            return kap_k, fg_k
    raise ValueError("no fg-kap or fg-fg interactions found")

def get_k_factor(l, type):
    """ return -1 if type not found """
    for x in l:
        if(x.type==type):
            return x.interaction_k_factor.value
    return -1
#    raise ValueError("type", type, "not found in list", l)


def round_way_up(val, bin_size = 1000):
    """
    rounds up to the floats array to the nearest multiple of
    an integral bin_size
    return value is an int
    """
    if not  isinstance(bin_size, numbers.Integral):
        raise TypeException("bin_size should be integral")
    ret = math.ceil(val / bin_size)  * bin_size
    return int(ret)

def print_list_to_file(FILEH, the_list, format="%d"):
    '''print list to filehandle FILEH'''
    first=True
    for t in the_list:
        if not first: FILEH.write(",")
        FILEH.write(format % t)
        first=False


def get_interaction_stats(interactions, is_a, is_b): #, is_sum_a=True):
    '''
    get the interaction statistics about interactions
    between types a and b, with statistics from point of view of a

    interatiocns - protobuf interaction statistics message
    is_a - boolean function to test if a type is a
    is_b - boolean function to test if a type is b
#    is_sum_a - if true, sum over all interactions that involve a, o/w avg

    return (pct_bnd,on_per_a_per_ns, off_per_a_per_ns)
    '''
    class istats:
        def __init__(self):
            self.pct_a = 0
            self.t=1E-10 # time
            self.on_per_a=0.0
            self.on_t=1E-10 # time for on
            self.off_per_a=0.0
            self.off_t=1E-10 # time for off
    itypes2stats={}
    for i in interactions:
        if is_a(i.type0) and is_b(i.type1):
            key=(i.type0, i.type1)
            if key not in itypes2stats:
                itypes2stats[key]=istats()
            s=itypes2stats[key]
            for op in i.order_params:
                # PCT BND
                d_t = op.misc_stats_period_ns
                s.pct_a = s.pct_a \
                    + d_t *  op.avg_fraction_bound_particles_i
                s.t = s.t + d_t
                # ON
                d_on_t = op.on_i_stats_period_ns
                s.on_per_a = s.on_per_a \
                    + d_on_t * op.avg_on_per_unbound_i_per_ns
                s.on_t = s.on_t + d_on_t
                # OFF
                d_off_t = op.off_stats_period_ns
                s.off_per_a = s.off_per_a \
                    + d_off_t * op.avg_off_per_bound_i_per_ns
                s.off_t = s.off_t + d_off_t
        if is_a(i.type1) and is_b(i.type0):
            key=(i.type0, i.type1)
            if key not in itypes2stats:
                itypes2stats[key]=istats()
            s=itypes2stats[key]
            for op in i.order_params:
                # PCT BND
                d_t = op.misc_stats_period_ns
                s.pct_a = s.pct_a + \
                    d_t * op.avg_fraction_bound_particles_ii
                s.t = s.t + d_t
                # ON
                d_on_t = op.on_ii_stats_period_ns
                s.on_per_a = s.on_per_a \
                    + d_on_t * op.avg_on_per_unbound_ii_per_ns
                s.on_t = s.on_t + d_on_t
                # OFF
                d_off_t = op.off_stats_period_ns
                s.off_per_a = s.off_per_a \
                    + d_off_t * op.avg_off_per_bound_ii_per_ns
                s.off_t = s.off_t + d_off_t
#    print t, on_t, off_t
    pct_a_per_t = sum(s.pct_a / s.t for s in itypes2stats.itervalues())
    on_per_a_per_t = sum(s.on_per_a / s.on_t for s in itypes2stats.itervalues())
    off_per_a_per_t = sum(s.off_per_a / s.off_t for s in itypes2stats.itervalues())
    return (pct_a_per_t, on_per_a_per_t, off_per_a_per_t)

def augment_results(pb_msg, results = {}):
    """
    read the protobuf message pb_msg and augment
    into the results dictionary

    returns true if message processed successfully
    """
    global KEY_CAPTIONS
    #  print a.is_valid()
    A = pb_msg.assignment
    S = pb_msg.statistics
    fg0_nbeads = A.fgs[0].number_of_beads.value
    fg0_R = A.fgs[0].radius.value
    for f in A.floaters:
        if is_kap(f.type):
            kap_R = f.radius.value
            kap_interactions = f.interactions.value
        if is_crap(f.type):
            crap_R = f.radius.value
    fgkap_k, fgfg_k = get_interaction_k( A.interactions )
    nonspecific_k = A.nonspecific_k.value
    for f in A.floaters:
        if is_kap(f.type):
            kap_k_factor = f.interaction_k_factor.value
    time_ns = round(S.bd_simulation_time_ns)
    time_step_fs = round(A.time_step)
    time_step_wave_factor = round(A.time_step_wave_factor.value)
    nup1_k_factor = get_k_factor( A.fgs, "Nup1_8copies" )
    rest_length_factor = A.fgs[0].rest_length_factor.value
    angular_D_factor = A.angular_D_factor.value
    work_unit = A.work_unit
    KEY_CAPTIONS = "kap_k_factor fgfg_k rest_length_factor kap_R nonspecific_k " \
                   "angular_D_factor time_step_fs time_step_wave_factor " \
                   + "nup1_k_factor fg0_R kap_interactions time_ns"
    key = tuple ( [ eval(k) for k in KEY_CAPTIONS.split() ] )
    if(not key in results):
        results[key] = {"n":0,
                        "kap_times":[],
                        "crap_times":[],
                        "small_crap_times":[],
                        "fg_length":[],
                        "kap_pct_bnd":[],
                        "crap_pct_bnd":[],
                        "representative_work_unit": work_unit
                    }
    print work_unit,
    results[key]["n"] = results[key]["n"] + 1
    if(A.slab_is_on.value):
        for a in S.floaters:
            if is_kap(a.type):
                print "KAP ",
                for t in a.transport_time_points_ns:
                    results[key]["kap_times"].append(t)
                    sys.stdout.write("%.1f," % t)
                sys.stdout.write(" ");
            if is_crap(a.type):
                print "INERT ",
                for t in a.transport_time_points_ns:
                    if(re.search("small", a.type)):
                        results[key]["small_crap_times"].append(t)
                        sys.stdout.write("%.1f*," % t)
                    else:
                        results[key]["crap_times"].append(t)
                        sys.stdout.write("%.1f," % t)
                sys.stdout.write(" ");
    print
    cyto_fgs = [
        "fg",
        "Nup100_8copies_chimera",
        "Nsp1_16copies_1",
        "Nup116_8copies_chimera"
    ];
    kap_pct, on_per_kap_per_ns, off_per_kap_per_ns= \
        get_interaction_stats(S.interactions,
                              is_kap,
                              is_fg)
                          #    lambda x: x in #cyto_fgs)
    crap_pct, on_per_crap_per_ns, off_per_crap_per_ns = \
        get_interaction_stats(S.interactions,
                              is_crap,
                              is_fg)
    results[key]["fg_length"].append(S.fgs[0].length)
    results[key]["kap_pct_bnd"].append(kap_pct)
    results[key]["crap_pct_bnd"].append(crap_pct)
    results[key]["kap_on"] = on_per_kap_per_ns
    results[key]["kap_off"] = off_per_kap_per_ns
    results[key]["crap_on"] = on_per_crap_per_ns
    results[key]["crap_off"] = off_per_crap_per_ns
    #        print key, results[key]


def print_results(results, file_name):
    FILE = open(file_name, "w")
    print >>FILE, KEY_CAPTIONS, "n transp_kaps transp_craps fg_length kap_pct_bnd crap_pct_bnd kap_on kap_off crap_on crap_off kap_transp_hist crap_transp_hist small_transp_hist representative_work_unit"
    for k, v in results.iteritems():
        for value in  k:
            print >>FILE, "%.2f" % (value),
        n=v["n"]
        print >>FILE, n,
        if "kap_times" in v:
            v["kap_times"].sort()
            print >>FILE, "%.2f" % (len(v["kap_times"]) * 1.0 / n),
        if "crap_times" in v:
            v["crap_times"].sort()
            print >>FILE, "%.2f" % (len(v["crap_times"]) * 1.0 / n),
        print >>FILE, "%.1f" % (sum(v["fg_length"]) / n),
        print >>FILE, "%.1f%%" % (100*sum(v["kap_pct_bnd"]) / n),
        print >>FILE, "%.1f%%" % (100*sum(v["crap_pct_bnd"]) / n),
        print >>FILE, "%.3f" % (v["kap_on"]),
        print >>FILE, "%.3f" % (v["kap_off"]),
        print >>FILE, "%.3f" % (v["crap_on"]),
        print >>FILE, "%.3f" % (v["crap_off"]),
        # print tranp histograms in bins
        bin_size = 5000
        if "kap_times" in v and "crap_times" in v:
            right_edge=round_way_up( max(v["kap_times"]+v["crap_times"] + v["small_crap_times"] +[bin_size]) , bin_size)
            bins=range(0,right_edge+bin_size,bin_size)
            h_kaps = numpy.histogram(v["kap_times"],bins)
            FILE.write(" ")
            print_list_to_file(FILE, h_kaps[0])
            h_craps=numpy.histogram(v["crap_times"],bins)
            FILE.write("  ")
            print_list_to_file(FILE, h_craps[0])
            h_small_craps=numpy.histogram(v["small_crap_times"],bins)
            FILE.write("  ")
            print_list_to_file(FILE, h_small_craps[0])
        FILE.write("  ")
        print >>FILE, v["representative_work_unit"],
        print >>FILE

def get_file_paths_generator(files_and_folders):
    """
    return a generator that iterates over:
    1) all the files that are not folders in files_and_folders
    2) all the inner *.pb files for every folder in files_and_folders
    """
    for path in files_and_folders:
        if(not os.path.exists(path)):
            print "WARNING: input path", path, "does not exist"
            continue
        if(os.path.isdir(path)):
            print "DIR PATH", path
            pb_paths = glob.glob(path + "/*.pb")
            for pb_path in pb_paths:
                print "PB_PATH", pb_path
                yield pb_path
        else:
            print "FILE_PATH", path
            yield path

def get_pb_messages_generator(files_and_folders):
    """
    return a generator that iterates over all pb output messages in:
    1) all the non-avro and avro files in files_and_folders
    2) all the inner *.pb for every folder in files_and_folders

    The generator returns a sequence of tuples (message, tag)
    with message the protobuf message, and the tag is the file name
    where the message is stored, with an internal id for avro files
    """
    o = Output()
    for file_path in get_file_paths_generator(files_and_folders):
        prefix,ext=os.path.splitext(file_path)
        if ext==".pb":
            try:
                FILE=open(file_path,"rb")
                o.ParseFromString(FILE.read())
                yield o, file_path
            except KeyboardInterrupt: exit(-1)
            except:
                print >> sys.stderr, 'Unexpected error: file %s' % file_path, sys.exc_info()
                continue
        else: # assume avro
            avro_reader=Avro2PBReader([file_path])
            i=0
            while(avro_reader.get_is_valid()):
                #  print a.is_valid()
                try:
                    i=i+1
                    s = avro_reader.read_next()
                    if(s == ""): break; # invalid output = the end
                    o.ParseFromString(s)
                    yield o, file_path + "." + str(i)
                except KeyboardInterrupt: exit(-1)
                except:
                    print >> sys.stderr, 'Unexpected error: file %s' % avro_reader.get_cur_file_name(), sys.exc_info()
                    print "CONTINUING"
                    continue






#################################
'''
Read pb or avro files in list of files and folders sys.arv[2:],
make stats on them, and write them to sys.argb[1]
'''
if __name__ != "__main__": exit()
output_file = sys.argv[1]
files_and_folders = sys.argv[2:]
print "Parsing %d files and folders" % len(files_and_folders)
if(os.path.exists(output_file)):
    if not query_yes_no("%s already exists - overwrite?" % output_file, None):
        exit(-1)
results = {}
n_total = 0
k_flush=10
for pb_message, tag in get_pb_messages_generator(files_and_folders):
    try:
        augment_results(pb_message, results)
        n_total = n_total + 1
        if(n_total % k_flush == 0):
            print "========== n = %d =========" % n_total
            print_results(results, output_file)
            pass
    except:
        type,msg,tb = sys.exc_info()
        print >> sys.stderr, 'Unexpected error in message tag %s - %s %s' \
            % (tag, type, msg)
        traceback.print_tb(tb)
# FINAL OUTPUT
print_results(results, output_file)
