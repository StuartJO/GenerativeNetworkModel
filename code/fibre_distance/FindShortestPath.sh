#!/bin/env bash
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH -t 2-0:0:0
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --export=ALL
#SBATCH --mem-per-cpu=32G
#SBATCH --account=kg98

% The directory where inputs/outputs from previous steps are stored, as well as where this scripts outputs will be stored
DATADIR=$1

% Directory where the code to run everything is
CODEDIR=$2

% The number of nodes to perform calculations of for each job
NODESSPERJOB=$3

% An optional fourth input gives the opportunity to redo/overwrite existing outputs if set to 1 
if [ -z "$4" ]; then
    REDO=$4
else
    REDO=0
fi

declare -i ID=$SLURM_ARRAY_TASK_ID

STARTNODE="$((1 + (${NODESSPERJOB} * ($ID - 1))))"
ENDNODE="$((${STARTPOINT} + ${NODESSPERJOB} - 1))"

module load matlab

# Unlike in other bits of code here, we don't just calculate the values in the upper triangle simply because I found there to be no computational benefit to doing this

if [ ${REDO} -eq 0 ]; then

if [ ! -f "${DATADIR}/DistanceMatrices/Output_${ID}.mat" ]; then

matlab -nodisplay -nosplash -r "addpath((genpath(${CODEDIR})); load('${DATADIR}/ReducedDirectMatrix.mat'); NODES = ${STARTNODE}:${ENDNODE}; NODES(NODES > Nodes)=[];for i = 1:length(POINTS); DIST(i,:) = distances(G,NODES(i)); end; save('${DATADIR}/DistanceMatrices/Output_${ID}.mat','DIST','NODES'); exit"

fi

else

matlab -nodisplay -nosplash -r "addpath((genpath(${CODEDIR})); load('${DATADIR}/ReducedDirectMatrix.mat'); NODES = ${STARTNODE}:${ENDNODE}; NODES(NODES > Nodes)=[];for i = 1:length(POINTS); DIST(i,:) = distances(G,NODES(i)); end; save('${DATADIR}/DistanceMatrices/Output_${ID}.mat','DIST','NODES'); exit"

fi
