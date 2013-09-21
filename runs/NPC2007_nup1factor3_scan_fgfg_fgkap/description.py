import glob
#configs=glob.glob("config*.pb")
#print configs
configs = ['config3_factor1_fg0.1.pb', 'config3_factor5_fg0.1.pb',
           'config3_factor1_fg1.0.pb', 'config3_factor5_fg1.0.pb',
           'config3_factor1_fg2.5.pb', 'config3_factor5_fg2.5.pb',
           'config3_factor1_fg3.5.pb', 'config3_factor5_fg3.5.pb',
           'config3_factor1_fg5.0.pb', 'config3_factor5_fg5.0.pb',
           'config3_factor1_fg10.0.pb', 'config3_factor5_fg10.0.pb']
binary="fg_simulation"
extra_arguments=["--short_init_factor=.1"]
distinct_work_units = 5
tasks_per_distinct = 2000
tasks=tasks_per_distinct * distinct_work_units* len(configs)
iterations=10
description=["""
10x2000s (=20 microsec) iterations of transport in coarse grained
NPC2007 model (FG motifs per bead = 2-3)
scanning mainly various time-step params and some others - mainly for run stability
"""]
