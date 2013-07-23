binary="fg_simulation"
extra_arguments=[]
tasks=30000
#tasks=2500
iterations=10
restart=["config_2.5_2.5",
         "config_3.0_3.5",
         "config_3.5_4.5",
         "config_4.5_3.0",
         "config_2.5_3.0",
         "config_3.0_4.0",
         "config_4.0_2.5",
         "config_4.5_3.5",
         "config_2.5_3.5",
         "config_3.0_4.5",
         "config_4.0_3.0",
         "config_4.5_4.0",
         "config_2.5_4.0",
         "config_3.5_2.5",
         "config_4.0_3.5",
         "config_4.5_4.5",
         "config_2.5_4.5",
         "config_3.5_3.0",
         "config_4.0_4.0",
         "config_3.0_2.5",
         "config_3.5_3.5",
         "config_4.0_4.5",
         "config_3.0_3.0",
         "config_3.5_4.0",
         "config_4.5_2.5"
         ]
description=["""
10x1000ns (=10 microsec) additional simulation time of transport through semi-real NPC
based on Nature 2007 model FG anchor points, starting from a simulation
that's been equilibrated for ~2.5 microsec, for a variety of k params
"""]
