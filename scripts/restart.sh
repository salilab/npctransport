#!/usr/bin/tcsh
set restart=$1
set k=$2
source ~/.aliases
npc_fg_simulation --restart $restart --conformations movie_$k.rmf --output out_${restart}_$k.pb --final_conformations final_${restart}_$k.rmf $3 $4 $5 $6 $7
