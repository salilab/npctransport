set cfg=$1
set k=$2
source ~/.aliases
npc_fg_simulation --configuration $cfg --conformations movie_$k.rmf --short_init_factor 0.05 --output out_$k.pb --final_conformations final_$k.rmf $3 $4 $5 $6 $7
