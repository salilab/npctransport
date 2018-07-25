#!/usr/bin/env python
from IMP.npctransport import *
import IMP.rmf
import RMF
import sys
f=open(sys.argv[1], "rb")
o= Output()
#config.SetTotalBytesLimit(50000000)

o.ParseFromString(f.read())
bch= RMF.BufferConstHandle(o.rmf_conformation);
RMF.write_buffer(bch, sys.argv[2] )
