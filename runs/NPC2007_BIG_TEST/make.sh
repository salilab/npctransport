#!/bin/bash
source ~/.aliases
echo "file fgfg rest_length_factor kap_R" > CONFIGS
n=1
for i in 0.01 0.1 0.25 1.0 2.5; do
    for j in 1.00 1.25 1.5;  do
        file=config$n.pb
        echo $file $i $j 25 >> CONFIGS
        $HOME/imp_git/fast/setup_environment.sh python make_cfg.py $file 25 $i $j;
        ((n++))
    done;
done
