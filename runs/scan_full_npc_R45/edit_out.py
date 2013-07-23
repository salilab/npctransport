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

def change_simulation_time(a, t_ns):
    """
    a - protobuf assignment object
    t - new time in ns
    """
    old_t_ns = a.simulation_time_ns
    time_step_ps = a.time_step # in picosecond = 10E-6 ns
    ps_in_ns = 1000000
    a.simulation_time_ns=t_ns
    a.number_of_frames = int( math.ceil(t_ns * ps_in_ns / time_step_ps) )


# fetch params
# Usage: <cmd> <outfile> <newoutfle>
outfile = sys.argv[1]
newoutfile = sys.argv[2]
new_diffusers_R = 45 # A
f=open(outfile, "rb")
output= Output()
output.ParseFromString(f.read())
f.close()
assignment = output.assignment
change_simulation_time(assignment, 1000)
for floater in assignment.floaters:
    floater.radius.value = new_diffusers_R
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
