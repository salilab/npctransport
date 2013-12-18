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
import collections
KEY_CAPTIONS = None



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


def dump_results(gop_table, n):
    OUT = open("OP.txt","w")
    is_first = True
    tags = []
    for t in sorted(gop_table.iterkeys()):
        if is_first:
            print >>OUT,"t",
            for tag in sorted(gop_table[t].iterkeys()):
                print >>OUT,tag,
                tags.append(tag)
            print >>OUT
            is_first=False
        print >>OUT, "%10.0f" % t,
        for tag in tags:
            print >>OUT, "%8.2f  " % ((gop_table[t][tag]+0.0)/n),
        print >>OUT

#################################
'''
Read pb or avro files in list of files and folders sys.arv[2:],
make stats on them, and write them to sys.argb[1]
'''
gop_table={}
if __name__ != "__main__": exit()
files_and_folders = sys.argv[1:]
print "Parsing %d files and folders" % len(files_and_folders)
n=0
for pb_message, tag in get_pb_messages_generator(files_and_folders):
    s = pb_message.statistics
    t_old = -1.0
    for gop in s.global_order_params:
        t= gop.time_ns
        if t==t_old:
            continue # skip to prevent double count
        if t not in gop_table:
            gop_table[t] = { "kap_interacting_fraction":0.0,
                             "kaps_channel":0,
                             "inerts_channel":0,
                             "fg_length":0.0,
                             "fg_rg":0.0}
        for i in range(12):
            i1 = i/3
            i2 = i%3
            tag = "zr_hists%d_%d" % (i1,i2)
            if not tag in gop_table[t]:
                gop_table[t][tag] = 0.0
            gop_table[t][tag] += gop.zr_hists[i1].ints[i2]
        t_old = t
    for fop in s.floaters[0].order_params:
        t = fop.time_ns
        assert t in gop_table
        gop_table[t]["kap_interacting_fraction"] += fop.interacting_fraction
        gop_table[t]["kaps_channel"] += fop.n_z1 + fop.n_z2
    for fop in s.floaters[1].order_params:
        t = fop.time_ns
        gop_table[t]["inerts_channel"] += fop.n_z1 + fop.n_z2
    for fg in s.fgs:
        for fgop in fg.order_params:
            t = fgop.time_ns
            assert t in gop_table
            gop_table[t]["fg_length"] += fgop.length / len(s.fgs)
            gop_table[t]["fg_rg"] += fgop.radius_of_gyration / len(s.fgs)
    n=n+1
    if ( n % 5 == 0):
        print "====== %d ======" % n
        dump_results(gop_table, n)
print
print "FINAL"
dump_results(gop_table, n)
