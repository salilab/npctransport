#!/bin/tcsh
source ~/.aliases
npc_restart $1 --conformations movie.R$2_$1.rmf --output output.R$2_$1.pb --statistics_level NONE> & LOG.R$2_$1 &