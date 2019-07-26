#!/usr/bin/env python
from IMP.npctransport import *
import sys
config= Configuration()
with open(sys.argv[1], "rb") as f:
config.ParseFromString(f.read())
f.close()
print(config)
