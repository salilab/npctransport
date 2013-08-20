#!/bin/bash
source ~/.aliases
for i in 0.01 0.1 0.25 1.0 2.5; do
    for j in 1.00 1.25 1.5;  do
        $HOME/imp_git/fast/setup_environment.sh python make_cfg.py NPC2007_R25_fgfg${i}_rlf${j}.cfg 25 $i $j;
    done;
done
