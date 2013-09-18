configs=["config2_f1.pb",
         "config3_f1.pb",
         "config2_f3.pb",
         "config3_f3.pb",
         "config2_f5.pb",
         "config3_f5.pb"
]
binary="fg_simulation"
extra_arguments=["--short_init_factor=.1"]
distinct_work_units = 120
tasks_per_distinct = 500
tasks=tasks_per_distinct * distinct_work_units* len(configs)
iterations=20
description=["""
20x1000s (=20 microsec) iterations of transport in coarse grained
NPC2007 model (FG motifs per bead = 2-3)
scanning mainly various time-step params and some others - mainly for run stability
"""]
