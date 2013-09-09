configs=["coarse2.pb",
         "coarse3.pb"]
binary="fg_simulation"
extra_arguments=["--short_init_factor=.1"]
tasks=120000*len(configs)
iterations=20
description=["""
20x500s (=15 microsec) iterations of transport in coarse grained
NPC2007 model (FG motifs per bead = 2-3)
scanning mainly various time-step params and some others - mainly for run stability
"""]
