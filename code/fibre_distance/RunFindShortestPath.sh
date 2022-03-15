#!/bin/env bash

% The directory where inputs/outputs from previous steps are stored, as well as where this scripts outputs will be stored
DATADIR=$1

% Directory where the code to run everything is
CODEDIR=$2

% The number of nodes to perform calculations of for each job
NODESSPERJOB=$3

% An optional input, set to 1 if ReducedDirectMatrix.mat has been made to skip making it again (otherwise it will be overridden)
if [ -z "$4" ]; then
    ReducedDirectMatrixMade=$4
else
    ReducedDirectMatrixMade=0
fi

module load matlab

if [ ${ReducedDirectMatrixMade} == 0 ]; then

rm ${DATADIR}/ShortestPathJobs.txt
rm ${DATADIR}/ReducedDirectMatrix.mat

matlab -nodisplay -nosplash -r "addpath((genpath(${CODEDIR})); G = CombineDirectConnectionResults('${DATADIR}'); Nodes = size(G.Nodes,1); dlmwrite('${DATADIR}/ShortestPathJobs.txt',ceil(Nodes / ${NODESSPERJOB})); save('${DATADIR}/ReducedDirectMatrix.mat','G','Nodes','-v7.3'); exit"

fi

NJOBS=$(sed -n "1p" /projects/kg98/stuarto/Fetal_brain/CRL_Fetal_Brain_Atlas_2017/NewOptimFibreLength/${BRAIN}ShortestPathJobs.txt)

echo "${NJOBS} jobs needed"

mkdir ${DATADIR}/FindShortestPathSlurmOutputs

mkdir ${DATADIR}/DistanceMatrices

SIMULTANEOUS_JOBS=240

sbatch --array=1-${NJOBS}%${SIMULTANEOUS_JOBS} --outputs=${DATADIR}/FindShortestPathOutputs/ ${CODEDIR}/FindShortestPath.sh ${DATADIR} ${NODESSPERJOB} 0