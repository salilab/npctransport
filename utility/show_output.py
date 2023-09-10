#!/usr/bin/env python
from IMP.npctransport import *
import sys
config= Output()
with open(sys.argv[1], "rb") as f:
    config.ParseFromString(f.read())
    print(config)
