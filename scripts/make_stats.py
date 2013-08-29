# #!/netapp/sali/barak/Runs/NPC_BD/Cyl2/tools/imppy.sh /usr/bin/python
from IMP.npctransport import *
import numpy
import math
import sys
import re
import numbers
import os

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
    raise ValueException("no fg-kap or fg-fg interactions found")


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


def get_interaction_stats(interactions, is_a, is_b):
    '''
    get the interaction statistics about interactions
    between types a and b, with statistics from point of view of a

    interatiocns - protobuf interaction statistics message
    is_a - boolean function to test if a type is a
    is_b - boolean function to test if a type is b

    return (pct_bnd,on_per_a_per_ns, off_per_a_per_ns)
    '''
    pct_a=0.0
    t=1E-10 # time
    on_per_a=0.0
    on_t=1E-10 # time for on
    off_per_a=0.0
    off_t=1E-10 # time for off
    for i in interactions:
        if is_a(i.type0) and is_b(i.type1):
            for op in i.order_params:
                # PCT BND
                d_t = op.misc_stats_period_ns
                pct_a = pct_a \
                    + d_t *  op.avg_fraction_bound_particles_i
                t = t + d_t
                # ON
                d_on_t = op.on_i_stats_period_ns
                on_per_a = on_per_a \
                    + d_on_t * op.avg_on_per_unbound_i_per_ns
                on_t = on_t + d_on_t
                # OFF
                d_off_t = op.off_stats_period_ns
                off_per_a = off_per_a \
                    + d_off_t * op.avg_off_per_bound_i_per_ns
                off_t = off_t + d_off_t
        if is_a(i.type1) and is_b(i.type0):
            for op in i.order_params:
                # PCT BND
                d_t = op.misc_stats_period_ns
                pct_a = pct_a + \
                    d_t * op.avg_fraction_bound_particles_ii
                t = t + d_t
                # ON
                d_on_t = op.on_ii_stats_period_ns
                on_per_a = on_per_a \
                    + d_on_t * op.avg_on_per_unbound_ii_per_ns
                on_t = on_t + d_on_t
                # OFF
                d_off_t = op.off_stats_period_ns
                off_per_a = off_per_a \
                    + d_off_t * op.avg_off_per_bound_ii_per_ns
                off_t = off_t + d_off_t
#    print t, on_t, off_t
    return (pct_a / t, on_per_a / on_t, off_per_a / off_t)


def augment_results(files, results = {}, max_entries = 1000):
    """
    read up to max_entries entries from avro reader and augment them
    into the results dictionary

    returns the number of entries actually read
    """
    global KEY_CAPTIONS
    o = Output()
    n_entries_read = 0
    for file_name in files:
        #  print a.is_valid()
        try:
            FILE=open(file_name,"rb")
            o.ParseFromString(FILE.read())
        except KeyboardInterrupt: exit(-1)
        except:
            print >> sys.stderr, 'Unexpected error: file %s' % file_name, sys.exc_info()
            continue
        A = o.assignment
        S = o.statistics
        try:
            fg_nbeads = A.fgs[0].number_of_beads.value
            kap_R = A.floaters[0].radius.value
            crap_R = A.floaters[1].radius.value
            fgkap_K, fgfg_K = get_interaction_k( A.interactions )
            nonspecific_K = A.nonspecific_k.value
            for f in A.floaters:
                if is_kap(f.type):
                    kap_K_factor = f.interaction_k_factor.value
            time_ns = round(S.bd_simulation_time_ns)
            time_step_fs = round(A.time_step)
        except:
            print >> sys.stderr, 'Unexpected error: file %s' % file_name, sys.exc_info()
            continue
        KEY_CAPTIONS = "kap_K_factor fgfg_K kap_R nonspecific_K time_ns time_step_fs"
#       KEY_CAPTIONS = "fg_nbeads kap_R crap_R fgkap_K fgfg_K nonspecific_K"
        key = tuple ( [ eval(k) for k in KEY_CAPTIONS.split() ] )
        if(not key in results):
            results[key] = {"n":0,
                            "kap_times":[],
                            "crap_times":[],
                            "fg_length":[],
                            "kap_pct_bnd":[],
                            "crap_pct_bnd":[]}
        print file_name,
        results[key]["n"] = results[key]["n"] + 1
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
                    results[key]["crap_times"].append(t)
                    sys.stdout.write("%.1f," % t)
                sys.stdout.write(" ");
        print
        cyto_fgs = [
            "Nup100_8copies_chimera",
            #            "Nsp1_16copies_1",
            #"Nup116_8copies_chimera"
            ];
        kap_pct, on_per_kap_per_ns, off_per_kap_per_ns= \
            get_interaction_stats(S.interactions,
                                  is_kap,
                                  lambda x: x in cyto_fgs)
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
        n_entries_read += 1
    return n_entries_read

def print_results(results, FILE=sys.stdout):
    print >>FILE, KEY_CAPTIONS, "n transp_kaps tranp_craps fg_length kap_pct_bnd crap_pct_bnd kap_transp_hist crap_transp_hist"
    for k, v in results.iteritems():
        for value in  k:
            print >>FILE, "%.2f" % (value),
        n=v["n"]
        print >>FILE, n,
        v["kap_times"].sort()
        v["crap_times"].sort()
        print >>FILE, "%.2f" % (len(v["kap_times"]) * 1.0 / n),
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
        right_edge=round_way_up( max(v["kap_times"]+v["crap_times"] +[bin_size]) , bin_size)
        bins=range(0,right_edge+bin_size,bin_size)
        h_kaps = numpy.histogram(v["kap_times"],bins)
        h_craps=numpy.histogram(v["crap_times"],bins)
        FILE.write(" ")
        print_list_to_file(FILE, h_kaps[0])
        FILE.write("  ")
        print_list_to_file(FILE, h_craps[0])
        print >>FILE

#################################
if __name__ != "__main__": exit()
output_file = sys.argv[1]
files=(sys.argv[2:])
if(os.path.exists(output_file)):
    if not query_yes_no("%s already exists - overwrite?" % output_file, None):
        exit(-1)
results = {}
n_total = 0
i=0
k=10
while(True):
    n_read=augment_results(files[i:i+k], results, max_entries = k)
    i=i+k
    n_total = n_total + n_read
    if(n_total == 0):
        print "NOTHING TO READ"
        break
#    print "========== n = %d =========" % n_total
    OUTPUT = open(output_file,"w")
    print_results(results, OUTPUT)
    del OUTPUT
    if(n_read == 0):
        break
