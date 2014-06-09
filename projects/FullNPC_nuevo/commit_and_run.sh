#!/bin/tcsh
if($#argv < 1) then
    echo need commit message
    exit
endif
exit

source ~/.aliases
git commit make_cfg.py commit_and_run.sh run.sh -m "$1"
/bin/tcsh -c run.sh
