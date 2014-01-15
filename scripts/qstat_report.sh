#!/bin/bash
for i in `qstat | awk '{n=n+1; if(n>2){print $1}}' | sort -u` ; do
    echo $i
    echo -n $i `qstat | grep $i | wc -l` " ";
    qstat -j $i | grep args;
done
