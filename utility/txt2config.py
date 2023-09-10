#!/usr/bin/env python
from IMP.npctransport import *
import google.protobuf.text_format
import sys
config= Configuration()
with open(sys.argv[1], "rb") as f:
    google.protobuf.text_format.Merge(f.read(),config)
with open(sys.argv[2], "wb") as fout:
    fout.write(config.SerializeToString())
