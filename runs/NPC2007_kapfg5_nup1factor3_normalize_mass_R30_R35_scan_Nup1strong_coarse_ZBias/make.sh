#!/bin/tcsh
set SET_IMP=~/imp_git/fast/setup_environment.sh
$SET_IMP python make_cfg.py config2_f1_R30.pb 30.0 2.0 1.0
$SET_IMP python make_cfg.py config3_f1_R30.pb 30.0 3.0 1.0
$SET_IMP python make_cfg.py config2_f3_R30.pb 30.0 2.0 3.0
$SET_IMP python make_cfg.py config3_f3_R30.pb 30.0 3.0 3.0
$SET_IMP python make_cfg.py config2_f5_R30.pb 30.0 2.0 5.0
$SET_IMP python make_cfg.py config3_f5_R30.pb 30.0 3.0 5.0
$SET_IMP python make_cfg.py config2_f1_R35.pb 35.0 2.0 1.0
$SET_IMP python make_cfg.py config3_f1_R35.pb 35.0 3.0 1.0
$SET_IMP python make_cfg.py config2_f3_R35.pb 35.0 2.0 3.0
$SET_IMP python make_cfg.py config3_f3_R35.pb 35.0 3.0 3.0
$SET_IMP python make_cfg.py config2_f5_R35.pb 35.0 2.0 5.0
$SET_IMP python make_cfg.py config3_f5_R35.pb 35.0 3.0 5.0
