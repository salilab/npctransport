#!/bin/csh
#$ -S /bin/csh
#$ -o /scrapp/barak/Logs/
#$ -cwd
#$ -j y
#$ -r y
#$ -N npc_default_restart
#$ -l arch=linux-x64,mem_free=0.75G
#$ -l h_rt=240:00:00
#$ -t 1-500

if($#argv < 1) then
    echo "Missing params running $0 $argv"
    echo
    echo "Usage: $0 restart_file [sim_time_factor] [out_folder_name]"
    echo "\n\t restart_file - file to start from"
    echo   "\t sim_time_factor - by how long to extend or shorten the simulation (default 1.0)"
    echo   "\t out_folder_name - note this is a relative location. By default,"
    echo   "\t                   it is the last relative part of the rmf_folder"
    exit -1
endif
echo "Running $0 $argv"

# Param setting
set IMP=/netapp/sali/barak/imp_git/fast/
#set NPC=/netapp/sali/barak/npc/fast/
set NPCBIN=$IMP/module_bin/npctransport/
set PYTHON=/netapp/sali/barak/MyPython/bin/python
set OUT=/scrapp/barak/
set CYL2=/netapp/sali/barak/Runs/NPC_BD/Cyl2/
set TMPDIR=/scratch
set MYTMP=`mktemp -d`
set seed=`od -An -N4 -td4 /dev/random`
echo SEED: $seed
set RESTART=`echo $1 | sed 's/\/$//g'`
if( ! -e $RESTART ) then
    echo $RESTART does not exist
    exit -1
endif
if($#argv >= 2) then
    set SIM_TIME_FACTOR=$2
else
    set SIM_TIME_FACTOR=1.0
endif
if($#argv >= 3) then
    set OUTFOLDER=$OUT/$3
else
    set OUTFOLDER=$OUT/${RESTART:t}_next
endif
set i=$SGE_TASK_ID
if(-d $RESTART) then
    set RESTARTFILE = $RESTART/out$i.pb
else
    set RESTARTFILE = $RESTART
endif
if(! -e $RESTARTFILE) then
    echo Restart file $RESTARTFILE does not exist
    exit -1
endif
echo "Restart file $RESTARTFILE"
echo "Output folder $OUTFOLDER"

if(! -e $OUT) mkdir $OUT
cd $OUT
if(! -e $OUTFOLDER) mkdir ${OUTFOLDER}
if(-e $OUTFOLDER/final$i.rmf) then
    echo Aborting: $OUTFOLDER/final$i.rmf exists
    exit -1
endif
if(! -e $OUTFOLDER/TIMESTAMP) then
    echo Job id: $JOB_ID > $OUTFOLDER/TIMESTAMP
endif
# Run:
cd $MYTMP
echo "Temporary run folder $MYTMP"
$IMP/setup_environment.sh $NPCBIN/fg_simulation --output $MYTMP/out$i.pb --conformations $MYTMP/movie$i.rmf --final_conformations $MYTMP/final$i.rmf --work_unit $i --random_seed $seed --restart $RESTARTFILE --short_sim_factor $SIM_TIME_FACTOR
set err=$status
if($err) then
    echo Error during run of fg_simulation - status code $err
    exit $err
endif
echo "Moving final output files from $MYTMP to $OUTFOLDER"
mv $MYTMP/* $OUTFOLDER
rmdir $MYTMP
