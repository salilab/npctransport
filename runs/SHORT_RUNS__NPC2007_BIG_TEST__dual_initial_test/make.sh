#!/bin/bash
source ~/.aliases
echo "file kaps_R coarse_factor" > CONFIGS
n=1
for i in 20 30; do # kaps_R
    for j in 1.0 2.0 3.0 ; do # coarse factor
        file=config$n.pb
        echo $file $i $j >> CONFIGS
        $HOME/imp_git/fast/setup_environment.sh python make_cfg.py $file $i $j
        ((n++))
    done;
done
