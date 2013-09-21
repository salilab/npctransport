#!/bin/tcsh
set SET_IMP=~/imp_git/fast/setup_environment.sh
foreach i (0.1 1.0 2.5 3.5 5.0 10.0)
    echo FG $i
    $SET_IMP python make_cfg.py config3_factor1_fg$i.pb 30.0 3.0 1.0 $i
    $SET_IMP python make_cfg.py config3_factor5_fg$i.pb 30.0 3.0 5.0 $i
end
