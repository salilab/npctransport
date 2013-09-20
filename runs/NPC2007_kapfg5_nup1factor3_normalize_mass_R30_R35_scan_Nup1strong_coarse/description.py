configs=["config2_f1_R30.pb",
         "config3_f1_R30.pb",
         "config2_f3_R30.pb",
         "config3_f3_R30.pb",
         "config2_f5_R30.pb",
         "config3_f5_R30.pb",
         "config2_f1_R35.pb",
         "config3_f1_R35.pb",
         "config2_f3_R35.pb",
         "config3_f3_R35.pb",
         "config2_f5_R35.pb",
         "config3_f5_R35.pb"
]
binary="fg_simulation"
extra_arguments=["--short_init_factor=.1"]
distinct_work_units = 4
tasks_per_distinct = 2000
tasks=tasks_per_distinct * distinct_work_units* len(configs)
iterations=20
description=["""
20x1000s (=20 microsec) iterations of transport in coarse grained
NPC2007 model (FG motifs per bead = 2-3)
scanning mainly various time-step params and some others - mainly for run stability
"""]
