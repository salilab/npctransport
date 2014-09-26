#!/usr/bin/env python
from IMP.npctransport import *
import sys
f=open(sys.argv[1], "rb")
output= Output()
output.ParseFromString(f.read())
f.close()
output.ClearField("rmf_conformation")
fo=open("NR" + sys.argv[1], "wb")
fo.write(output.SerializeToString())
fo.close
