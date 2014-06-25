#!/usr/bin/python
from IMP.npctransport import *
import sys
f=open(sys.argv[1], "rb")
config= Output()
config.ParseFromString(f.read())
print config.statistics.floaters[1]
print config.statistics.floaters[2]
