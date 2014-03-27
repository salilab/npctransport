#!/usr/bin/python
from IMP.npctransport import *
import sys
import math

outfile = sys.argv[1]
newoutfile = sys.argv[2]
f=open(outfile, "rb")
output= Output()
output.ParseFromString(f.read())
f.close()
a = output.assignment
for f in a.floaters:
    f.k_z_bias.value=-1.0
    f.k_z_bias_fraction.value=1.0
# dump to new file
nf=open(newoutfile, "wb")
nf.write(output.SerializeToString())
print output.assignment
