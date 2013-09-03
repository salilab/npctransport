#!/usr/bin/tcsh
set restart=$1
source ~/.aliases
npc_fg_simulation --restart $restart --conformations movie_${restart:r}.rmf --output out_${restart:r}.pb --final_conformations final_${restart:r}.rmf $2 $3 $4 $5 $6 $7 >& LOG.${restart:r}
