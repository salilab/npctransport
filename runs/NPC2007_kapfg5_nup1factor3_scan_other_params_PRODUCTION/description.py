configs=["coarse2.pb",
         "coarse3.pb"]
binary="fg_simulation"
extra_arguments=["--short_init_factor=.1"]
#tasks=120000*len(configs)
distinct_work_units = 120
tasks_per_distinct = 1000
tasks=tasks_per_distinct * distinct_work_units* len(configs)
iterations=40
description=["""
40x375s (=15 microsec) iterations of transport in coarse grained
NPC2007 model (FG motifs per bead = 2-3)
scanning mainly various time-step params and some others - mainly for run stability
"""]
