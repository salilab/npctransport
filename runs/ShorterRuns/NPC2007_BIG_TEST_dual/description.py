configs=[]
for i in range(6):
    configs.append("conifg%d.pb" % (i+1))
print configs
binary="fg_simulation"
#extra_arguments=["--short_init_factor=.1"]
#tasks=4147200*len(configs) # for production
tasks=1000*len(configs) #for testing
iterations=2
description=["""
massive scan of different params
"""]
