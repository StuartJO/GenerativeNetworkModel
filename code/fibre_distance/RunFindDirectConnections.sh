#!/bin/env bash

# This is just a wrapper script for FindDirectConnections.sh

% The directory where inputs/outputs from previous steps are stored, as well as where this scripts outputs will be stored
DATADIR=$1

% Directory where the code to run everything is
CODEDIR=$2

% The number of vertices/points to perform calculations of for each job
POINTSPERJOB=$3

NJOBS=$(sed -n "1p" ${DATADIR}/ReducedJobs.txt)

echo "${NJOBS} jobs needed"

% 240 CPUS is the maximum I had access to!
SIMULTANEOUS_JOBS=240

mkdir ${DATADIR}/FindDirectConnectionsSlurmOutputs

mkdir ${DATADIR}/DirectConnections

sbatch --array=1-${NJOBS}%${SIMULTANEOUS_JOBS} --output=${DATADIR}/FindDirectConnectionsSlurmOutputs ${CODEDIR}/FindDirectConnections.sh ${DATADIR} ${CODEDIR} ${POINTSPERJOB}
