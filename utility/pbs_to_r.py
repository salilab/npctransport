#!/usr/bin/env python
import sys
import IMP.npctransport

"""Take a list of protobufs from assignments and statistics, assuming they are
named like foo.statistics.pb and foo.assignments.pb, extract the listed fields
and print the results in a file suitable for input into r"""

def fix(s):
  stage1= s.replace(".", "_").replace("(", "").replace(")", "")
  if stage1.endswith("_value"):
    return stage1[:-len("_value")]
  else:
    return stage1

assignments_fields=["nonspecific_range.value", "interaction_spring_constant.value"]
statistics_fields=["energy_per_particle"]
data={}

for f in sys.argv[1:]:
  # silly assumption about no . in dir names
  pref=f[f.rfind("/")+1:f.find(".")]
  if pref not in data.keys():
    data[pref]={}
  #print f
  if f.endswith(".assignments.pb"):
    pb= IMP.npctransport.Assignment()
    fields=assignments_fields
  elif f.endswith(".statistics.pb"):
    pb= IMP.npctransport.Statistics()
    fields=statistics_fields
  pb.ParseFromString(open(f, "rb").read())
  for a in fields:
    call="pb."+a
    vs=eval(call)
    #print call,vs
    v=float(vs)
    data[pref][fix(a)]=v

allf= [fix(x) for x in assignments_fields+statistics_fields]
print "work_unit",
for a in allf:
  print a,
print
for d in data.keys():
  print d,
  for a in allf:
    print data[d][a],
  print
