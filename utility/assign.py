#!/usr/bin/env python

import IMP.npctransport
import sys
try:
    import argparse
except ImportError:
    import IMP._compat_argparse as argparse

p = argparse.ArgumentParser()
p.add_argument("-w", "--work_unit", dest="work_unit", type=int,
               help="The work unit")
p.add_argument("configuration", help="configuration file name")
p.add_argument("assignment", help="assignment file name")

args = p.parse_args()
if not args.work_unit:
    p.print_help()
    sys.exit(1)
config= IMP.npctransport.Configuration()
IMP.npctransport.assign_ranges(args.configuration, args.assignment,
                               args.work_unit, True)
