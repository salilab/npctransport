#!/usr/bin/env python
from IMP.npctransport import *
import google.protobuf.text_format
import sys
f=open(sys.argv[1], "rb")
config= Configuration()
google.protobuf.text_format.Merge(f.read(),config)
f.close()
#print config
# dump to file
fout=open(sys.argv[2], "wb")
fout.write(config.SerializeToString())
fout.close()
