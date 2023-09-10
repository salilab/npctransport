#!/usr/bin/env python
from IMP.npctransport import *
import sys
with open(sys.argv[1], "rb") as f:
	output= Output()
	output.ParseFromString(f.read())
print(output.assignment)
