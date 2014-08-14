#!/bin/tcsh
source ~/.aliases
set output=output.R$2_$1.pb
if ( -e $output ) then
  echo Output file $output already exists - delete to restart
  exit -1
endif
npc_restart $1 --conformations movie.R$2_$1.rmf --output output.R$2_$1.pb --statistics_level NONE --short_sim_factor 4.0 > & LOG.R$2_$1 &
