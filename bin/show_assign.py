#!/usr/bin/env python
from IMP.npctransport import *
import sys
f=open(sys.argv[1], "rb")
output= Output()
output.ParseFromString(f.read())
print output.assignment
