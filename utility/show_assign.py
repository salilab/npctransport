#!/usr/bin/env python
from IMP.npctransport import *
import sys
output= Output()
with open(sys.argv[1], "rb") as f:
    output.ParseFromString(f.read())
print(output.assignment)

