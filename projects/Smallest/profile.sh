#!/bin/tcsh
if($#argv < 1) then
    echo missing output restart file name
    exit -1
endif
source ~/.aliases
npc_restart $1 --cpu_profile --statistics_level NONE --short_sim_factor 0.0075
