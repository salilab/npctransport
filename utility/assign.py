#!/usr/bin/env python

import IMP.npctransport
import sys
import optparse

parser = optparse.OptionParser(usage="usage: %prog configuration.pb assignment.pb")
parser.add_option("-w", "--work_unit", dest="work_unit",
                  help="The work unit")
(options, args) = parser.parse_args()
if len(args) != 2 or not options.work_unit:
        parser.print_help()
        exit(1)
config= IMP.npctransport.Configuration()
IMP.npctransport.assign_ranges(args[0], args[1], int(options.work_unit), True)
