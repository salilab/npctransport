#!/usr/bin/python
from IMP.npctransport import *
import sys
f=open(sys.argv[1], "rb")
output= Output()
output.ParseFromString(f.read())
assign = output.assignment
print "k", assign.interactions[0].interaction_k
fg = assign.fgs[0]
print "tamd params", fg.is_tamd, \
    fg.tamd_T_factor_coeff.value, \
    fg.tamd_T_factor_base.value, \
    fg.tamd_T_factor_coeff.value, \
    fg.tamd_T_factor_base.value, \
    fg.tamd_K.value
print output.statistics.floaters[0]
print output.statistics.floaters[1]
print output.statistics.floaters[2]
