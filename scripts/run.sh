#!/usr/bin/tcsh
set OUTFOLDER=Tmp
set cfg=$1
set ext=$2
set sif=$3
if($#argv < 3) then
    set sif=0.05
endif
source ~/.aliases
if(! -e $OUTFOLDER) then
    mkdir $OUTFOLDER
endif
if(! -d $OUTFOLDER) then
    echo $OUTFOLDER already exists as file, choose a different name
    exit -1
endif
echo npc_fg_simulation --configuration $cfg --conformations Tmp/movie_$ext.rmf --short_init_factor $sif \
    --output Tmp/out_$ext.pb --final_conformations Tmp/final_$ext.rmf $4 $5 $6 $7
npc_fg_simulation --configuration $cfg --conformations Tmp/movie_$ext.rmf --short_init_factor $sif \
    --output Tmp/out_$ext.pb --final_conformations Tmp/final_$ext.rmf $4 $5 $6 $7
