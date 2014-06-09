#!/bin/tcsh
source ~/.aliases
if ( `git status ./make_cfg.py run2.sh -s | wc -w` > 0 ) then
    echo commit make_cfg.py and run.sh please
    exit -1
endif
set k=`git log --pretty=format:"%h" -n 1`
echo GIT revision $k;
imppy make_cfg.py test_$k.cfg;
npc_fg_simulation test_$k.cfg --cylinder_nlayers 3 --conformations movie_$k.rmf --short_sim_factor 0.1 > & LOG_$k &
