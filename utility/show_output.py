#!/usr/bin/env python
from IMP.npctransport import *
import sys
config= Output()
#config.SetTotalBytesLimit(50000000)
with open(sys.argv[1], "rb") as f:
    config.ParseFromString(f.read())
    #line=f.read(100)
    #while line<>"":
    #    config.MergeFromString(line)
    #    line= f.read(100)
    print(config)
