#!/usr/bin/tcsh
set cfg=$1
set k=$2
set sif=$3
if($#argv < 3) then
    set sif=0.05
endif
source ~/.aliases
echo npc_fg_simulation --configuration $cfg --conformations movie_$k.rmf --short_init_factor $sif --output out_$k.pb --final_conformations final_$k.rmf $3 $4 $5 $6 $7
 npc_fg_simulation --configuration $cfg --conformations movie_$k.rmf --short_init_factor $sif --output out_$k.pb --final_conformations final_$k.rmf $4 $5 $6 $7
