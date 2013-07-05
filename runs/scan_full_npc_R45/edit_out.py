#!/usr/bin/python
from IMP.npctransport import *
import sys
import math
import re
from numpy import arange

def is_fg_kap(i):
    """ checks if interaction i is between fg and kap """
    return ( re.match("fg", i.type0) and re.match("kap", i.type1) ) \
        or ( re.match("kap", i.type0) and re.match("fg", i.type1) )

def is_fg_fg(i):
    """ checks if interaction i is between fg and fg """
    return ( re.match("fg", i.type0) and re.match("fg", i.type1) ) \

# fetch params
# Usage: <cmd> <outfile> <newoutfle>
outfile = sys.argv[1]
newoutfile = sys.argv[2]
f=open(outfile, "rb")
output= Output()
output.ParseFromString(f.read())
f.close()
assignment = output.assignment
assignment.simulation_time_ns=1000
for floater in assignment.floaters:
    floater.radius.value = 45
for k1 in arange(2.5,4.51,0.5):
    for k2 in arange(2.5,4.51,0.5):
        for i in assignment.interactions:
            if (is_fg_kap(i)) :
                i.interaction_k.value = k1
            if (is_fg_fg(i)) :
                i.interaction_k.value = k2
        # dump to new file
        newoutfile_kk = "%s_%.1f_%.1f" % (newoutfile, k1, k2)
        nf=open(newoutfile_kk, "wb")
        nf.write(output.SerializeToString())
