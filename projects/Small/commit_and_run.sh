#!/bin/tcsh

if ( `git status ./make_cfg.py run2.sh -s | wc -w` > 0 ) then
    if($#argv < 1) then
        echo need commit message
        exit
    endif
    echo Committing with message \'$1\'
    git status -uno
    git add ./make_cfg.py ./commit_and_run.sh ./run.sh
    git commit  make_cfg.py commit_and_run.sh run.sh -m "$1"
else
    echo Nothing to commit - running
endif
/bin/tcsh -c ./run.sh
