#!/bin/tcsh
set SET_IMP=~/imp_git/fast/setup_environment.sh
$SET_IMP python make_cfg.py coarse2.pb 30.0 2.0
$SET_IMP python make_cfg.py coarse3.pb 30.0 3.0
