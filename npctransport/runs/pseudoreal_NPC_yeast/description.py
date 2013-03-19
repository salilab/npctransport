binary="fg_simulation"
extra_arguments=["--cylinder_nlayers", "5"]
tasks=400000 # for 20x20 steps x 1000 jobs
#tasks=2500
iterations=8
description=["""
5x2000ns (=10microsec) simulations of transport aimed for scanning:
* non-specific interaction
* specific interaction of fg-fg and fg-TF
for a semi-realistic NPC (thought of right dimenstions)
"""]
