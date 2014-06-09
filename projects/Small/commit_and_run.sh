#!/bin/tcsh
if($#argv < 1) then
    echo need commit message
    exit
endif


echo Committing with message \'$1\'
git add  make_cfg.py commit_and_run.sh run.sh
git status -uno
git commit -m "$1"
/bin/tcsh -c run.sh
