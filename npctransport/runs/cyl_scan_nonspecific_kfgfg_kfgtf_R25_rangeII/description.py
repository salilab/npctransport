binary="fg_simulation"
extra_arguments=["--cylinder_nlayers", "4"]
tasks=300000 # for 3x50x20 steps x 150 jobs
#tasks=2500
iterations=8
description=["""
8x2500ns (=20ms) simulations of transport aimed for scanning:
* non-specific interaction
* specific interaction of fg-fg and fg-TF
the range used is lower fg-fg interaction strength,
and denser than in previous tests
"""]
