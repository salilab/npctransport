#!/usr/bin/env python
from IMP.npctransport import *
import sys
f=open(sys.argv[1], "rb")
config= Configuration()
config.ParseFromString(f.read())
f.close()
print config
