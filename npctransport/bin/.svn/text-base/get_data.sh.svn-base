#!/bin/zsh
name=$1
iteration=$2
dest=~/results/$1
mkdir -p $dest
mkdir -p $dest/final
mkdir -p $dest/output
versname="$name"_"$iteration"
#f cp -parallel_copy=15 $CHOME$versname/output/\*.out $dest/output
cat `ls $dest/output/*.out | head -n 1` | head -n 1 > $dest/results.csv
cat $dest/output/*.out | grep -v "^#" >> $dest/results.csv
f cp -parallel_copy=15 $CHOME$versname/output/\*.final.pym $dest/final/