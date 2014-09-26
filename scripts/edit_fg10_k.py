#!/usr/bin/env python
from IMP.npctransport import *
import sys
import math
#import read_nups

# fetch params
# changes inoutfile in place with k=k*new_k_factor for fg10 (Nup1)
# Usage: <cmd> <inoutfile> <new_k_factor>
inoutfile = sys.argv[1]
new_k_factor = float(sys.argv[2])
f=open(inoutfile, "rb")
output= Output()
output.ParseFromString(f.read())
f.close()
assignment = output.assignment
for interaction in assignment.interactions:
    if(interaction.type0 == "fg10" and interaction.type1 == "kap"):
        old_k = interaction.interaction_k.value
        interaction.interaction_k.value *= new_k_factor
        print "Switching from old ", old_k, " to ", interaction.interaction_k.value
        break
# dump
nf=open(inoutfile, "wb")
nf.write(output.SerializeToString())
