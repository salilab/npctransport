#!/usr/bin/env python

import RMF
import IMP.npctransport
import sys

pb = IMP.npctransport.Output()
pb.ParseFromString(open(sys.argv[1], "rb").read())
bch = RMF.BufferConstHandle(pb.rmf_conformation)
RMF.write_buffer(bch, sys.argv[2])


# test it
bch = RMF.read_buffer(sys.argv[2])
fh = RMF.open_rmf_buffer_read_only(bch)
