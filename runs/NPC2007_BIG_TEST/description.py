configs=["conifg1.pb",
         "conifg2.pb",
         "conifg3.pb",
         "conifg4.pb",
         "conifg5.pb",
         "conifg6.pb",
         "conifg7.pb",
         "conifg8.pb",
         "conifg9.pb",
         "conifg10.pb",
         "conifg11.pb",
         "conifg12.pb",
         "conifg13.pb",
         "conifg14.pb",
         "conifg15.pb"]
binary="fg_simulation"
extra_arguments=["--short_init_factor=.1"]
tasks=22500*len(configs)
#300000
#tasks=2500
iterations=40
description=["""
40x500ns (=20 microsec) iterations of transport in coarse grained
NPC2007 model (FG motifs per bead = 3)
scanning:
1) lotta stuff
"""]
