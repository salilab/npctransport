#!/usr/bin/env python
from IMP.npctransport import *
import sys
f=open(sys.argv[1], "rb")
config= Output()
config.ParseFromString(f.read())
print(config.statistics)
