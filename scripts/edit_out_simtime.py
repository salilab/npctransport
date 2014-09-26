#!/usr/bin/env python
from IMP.npctransport import *
import sys
import math
#import read_nups

inoutfile = sys.argv[1]
new_t_ns = float(sys.argv[2])
f=open(inoutfile, "rb")
output= Output()
output.ParseFromString(f.read())
f.close()
assignment = output.assignment
old_t_ns = assignment.simulation_time_ns
old_nf = assignment.number_of_frames
assignment.simulation_time_ns = new_t_ns
new2old = (new_t_ns + 0.0) / old_t_ns
assignment.number_of_frames = int( math.ceil ( new2old * old_nf ) )
print "new # of frames", assignment.number_of_frames
# dump
nf=open(inoutfile, "wb")
nf.write(output.SerializeToString())
#print output.assignment.floaters
