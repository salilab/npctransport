#!/bin/tcsh
source ~/.aliases
if ( `git status ./make_cfg.py run2.sh -s | wc -w` > 0 ) then
    echo commit make_cfg.py and run.sh please
    exit -1
endif
set R=`git log --pretty=format:"%h" -n 1`
echo GIT revision $R;
set k=${R}$1
imppy make_cfg.py test_$k.cfg;
set init_only="--short_sim_factor 0.00001"
npc_fg_simulation test_$k.cfg --conformations movie_$k.rmf --short_init_factor 0.1 $init_only  --output output_$k.pb --statistics_level NONE > & LOG_$k &
