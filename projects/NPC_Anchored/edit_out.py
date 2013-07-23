#!/usr/bin/python
from IMP.npctransport import *
import sys
import math
import read_nups

# fetch params
# Usage: <cmd> <outfile> <newoutfle>
outfile = sys.argv[1]
newoutfile = sys.argv[2]
new_time_ns = float(sys.argv[3])
new_dump_ns = float(sys.argv[4])
f=open(outfile, "rb")
output= Output()
output.ParseFromString(f.read())
f.close()
assignment = output.assignment
old_time_ns = assignment.simulation_time_ns
old_dump_ns = assignment.dump_interval_ns
assignment.simulation_time_ns = new_time_ns
assignment.dump_interval_ns = new_dump_ns
assignment.number_of_frames = \
    int((new_time_ns / old_time_ns) * assignment.number_of_frames)
assignment.dump_interval_frames = \
    int((new_dump_ns / old_dump_ns) * assignment.dump_interval_frames)
# dump to new file
nf=open(newoutfile, "wb")
nf.write(output.SerializeToString())
print output.assignment
