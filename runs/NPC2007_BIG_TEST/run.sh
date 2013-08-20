set cfg=$1
set k=$2
source ~/.aliases
set imp="$HOME/imp_git/fast/setup_environment.sh"
set npc_fg_simulation="$imp $HOME/imp_git/fast/module_bin/npctransport/fg_simulation"
$npc_fg_simulation --configuration $cfg --conformations movie_$k.rmf --short_init_factor 0.05 --output out_$k.pb --final_conformations final_$k.rmf $3 $4 $5 $6 $7
