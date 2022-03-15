#!/bin/env bash
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH -t 3-0:0:0
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --export=ALL
#SBATCH --mem-per-cpu=8G
#SBATCH --account=kg98

# This script is implicitly intended to work via sumbission via RunFindDirectConnections.sh

% The directory where inputs/outputs from previous steps are stored, as well as where this scripts outputs will be stored
DATADIR=$1

% Directory where the code to run everything is
CODEDIR=$2

% The number of vertices/points to perform calculations of for each job
POINTSPERJOB=$3

% An optional fourth input gives the opportunity to redo/overwrite existing outputs if set to 1 
if [ -z "$4" ]; then
    REDO=$4
else
    REDO=0
fi

declare -i ID=$SLURM_ARRAY_TASK_ID

STARTPOINT="$((1 + (${POINTSPERJOB} * ($ID - 1))))"
ENDPOINT="$((${STARTPOINT} + ${POINTSPERJOB} - 1))"

module load matlab

if [ ${REDO} -eq 0 ]; then

if [ ! -f "${DATADIR}/DirectConnections/Output_${ID}.mat" ]; then

matlab -nodisplay -nosplash -r "addpath(genpath('${CODEDIR}')); load('${DATADIR}/Reduced.mat'); source_inds = ${STARTPOINT}:${ENDPOINT}; [DIRECT,source_inds] = RunDirectConnectionSurface(surface,target_points,source_inds);save('${DATADIR}/DirectConnections/Output_${ID}.mat','source_inds','DIRECT'); exit"

fi

else

matlab -nodisplay -nosplash -r "addpath(genpath('${CODEDIR}')); load('${DATADIR}/Reduced.mat'); source_inds = ${STARTPOINT}:${ENDPOINT}; [DIRECT,source_inds] = RunDirectConnectionSurface(surface,target_points,source_inds);save('${DATADIR}/DirectConnections/Output_${ID}.mat','source_inds','DIRECT'); exit"

fi