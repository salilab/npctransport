#!/bin/csh
#$ -S /bin/csh
#$ -o /scrapp/barak/Logs/
#$ -cwd
#$ -j y
#$ -r y
#$ -N npc_default
#$ -l arch=linux-x64,mem_free=0.75G
#$ -l h_rt=240:00:00
#$ -t 1-500

if($#argv < 1) then
    echo "Usage: $0 <cfg_file> <sim_time_factor> [out_folder_name] [work_unit]"
    echo   "\t sim_time_factor - by how long to extend or shorten the simulation (default 1.0)"
    echo "\n\t out_folder_name - either a relative or absolute location"
    echo "\n\t work_unit - the work_unit to send [the default is the task numeer]
    exit -1
endif
echo Running \"$0 $argv\"

# Param setting
set IMP=/netapp/sali/barak/imp_git/fast/
#set NPC=/netapp/sali/barak/npc/fast/
set NPCBIN=$IMP/module_bin/npctransport/
set PYTHON=/netapp/sali/barak/MyPython/bin/python
set OUT=/scrapp/barak/
set TMPDIR=/scratch
set MYTMP=`mktemp -d`
set seed=`od -An -N4 -td4 /dev/random`
echo SEED: $seed
set cfg_full=$1
set cfg_file=$cfg_full:t
set cfg_id=$cfg_file:r
if($#argv >= 2) then
    set SIM_TIME_FACTOR=$2
else
    set SIM_TIME_FACTOR=1.0
endif
if($#argv >= 3) then
    set OUTFOLDER=$OUT/$3
else
    set OUTFOLDER=$OUT/${cfg_id}
endif
if($#argv >= 4) then
    set WORK_UNIT=$4
else
    set WORK_UNIT=$i
endif
set i=$SGE_TASK_ID
echo "Cfg file $cfg_file  ;  Work id $i"
echo "Output folder $OUTFOLDER"

# Outfolder:
if(! -e $OUT) mkdir $OUT
cd $OUT
if(! -e $OUTFOLDER) mkdir ${OUTFOLDER}
if(-e $OUTFOLDER/out$i.pb) then
    echo Aborting: $OUTFOLDER/out$i.pb exists
    exit -1
endif
if(! -e $OUTFOLDER/TIMESTAMP) then
    echo Job id: $JOB_ID > $OUTFOLDER/TIMESTAMP
endif

# Config:
set NEWCFG = $OUTFOLDER/$cfg_file
set NEWCFG_txt = $OUTFOLDER/$cfg_id.txt
cp $cfg_full $NEWCFG
$IMP/setup_environment.sh python $NPCBIN/show_config.py $cfg_full > $NEWCFG_txt

# Run:
cd $MYTMP
echo "Temporary run folder $MYTMP"
$IMP/setup_environment.sh $NPCBIN/fg_simulation --configuration $OUTFOLDER/$cfg_file --output $MYTMP/out$i.pb --conformations $MYTMP/movie$i.rmf --final_conformations $MYTMP/final$i.rmf --work_unit $WORK_UNIT --random_seed $seed --short_init_factor 0.2 --short_sim_factor $SIM_TIME_FACTOR
set err=$status
if($err) then
    echo Error during run of fg_simulation - status code $err
    exit $err
endif
echo "Moving final output files from $MYTMP to $OUTFOLDER"
mv $MYTMP/* $OUTFOLDER
rmdir $MYTMP
