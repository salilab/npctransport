#!/bin/tcsh
set SET_IMP=~/imp_git/fast/setup_environment.sh
$SET_IMP python make_cfg.py config0.5.pb 30.0 0.5
$SET_IMP python make_cfg.py config.pb 30.0 1.0
$SET_IMP python make_cfg.py config2.pb 30.0 2.0
$SET_IMP python make_cfg.py config3.pb 30.0 3.0
