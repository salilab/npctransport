#!/usr/bin/python
from IMP.npctransport import *
import google.protobuf.text_format
import sys
f=open(sys.argv[1], "rb")
o= Output()
google.protobuf.text_format.Merge(f.read(),o)
f.close()
# dump to file
fout=open(sys.argv[2], "wb")
fout.write(o.SerializeToString())
fout.close()
