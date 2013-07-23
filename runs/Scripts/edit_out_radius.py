#!/usr/bin/python
from IMP.npctransport import *
import sys
import math
#import read_nups

# fetch params
# changes inoutfile in place with new floaters radii
# Usage: <cmd> <inoutfile> <newoutfle>
inoutfile = sys.argv[1]
new_radius = float(sys.argv[2])
f=open(inoutfile, "rb")
output= Output()
output.ParseFromString(f.read())
f.close()
assignment = output.assignment
for floater in assignment.floaters:
    print floater
    old_R = floater.radius.value
    floater.radius.value = new_radius
# dump
nf=open(inoutfile, "wb")
nf.write(output.SerializeToString())
print output.assignment.floaters
