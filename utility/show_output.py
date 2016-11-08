#!/usr/bin/env python
from IMP.npctransport import *
import sys
f=open(sys.argv[1], "rb")
config= Output()
#config.SetTotalBytesLimit(50000000)

config.ParseFromString(f.read())
#line=f.read(100)
#while line<>"":
#    config.MergeFromString(line)
#    line= f.read(100)
print config
