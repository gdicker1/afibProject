#!/bin/sh
CMADIR="/gscratch/gdicker1/CMAES_runner"
RUNDIR="$CMADIR/Text_Files/homogeneous/$(date +%s)"
mkdir -p $RUNDIR
cd $RUNDIR
### load relevant modules
### module load singularity/2.3.2
### start the job
echo "Starting job at $(date)"
### Run the python script with the afib.simg container specifying
###  the directory of the batchtool, directory of the reference tissue
###  and the indices of the APD, Resistance, and Parameter files to modify
singularity exec -B $CMADIR $CMADIR/afib.simg  \
python $CMADIR/homog_runner.py $CMADIR/Batchtool/ \
$CMADIR/Tissues/ 10 19
echo "Job finished at $(date)"
