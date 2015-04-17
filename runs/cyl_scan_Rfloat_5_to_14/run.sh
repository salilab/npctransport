#!/bin/tcsh
source ~/.aliases
set IMP=$1
set IMPNPC=$IMP/module_bin/npctransport/
set i=$2
set seed=$3
$IMP/setup_environment.sh $IMPNPC/fg_simulation --configuration config.pb --cylinder_nlayers 4 --output output$i.pb --conformations movie$i.rmf --work_unit $i --short_init_factor 0.5 --short_sim_factor 1.0 --random_seed $seed
