#!/bin/env bash

# This script is designed to a) compute a reduce surface b) find points on the interior to the surface perpendicular to each vertex of the reduced surface c) submit jobs to find all posssible direct connections in the reduced surface

#for i in {37..38}; do ./INITALISE_SURFACE.sh $i 100 1; done

BRAIN=$1

POINTSPERJOB=$2
#POINTSPERJOB=100
REDO=$3

HEMI="l"
REGTYPE="sulc"
MSMFOLDER="MSM_outputs_smoothed_upsampled_hocr"
NAME=""
FACTOR=0.15;

SURFACE_DIR="/projects/kg98/stuarto/Fetal_brain/CRL_Fetal_Brain_Atlas_2017"
GIFTI_DIR="/projects/kg98/stuarto/Fetal_brain/CRL_Fetal_Brain_Atlas_2017/gifti-1.6/" 

module load matlab

if [ ${REDO} -eq 0 ]; then

rm -fv /projects/kg98/stuarto/Fetal_brain/CRL_Fetal_Brain_Atlas_2017/NewOptimFibreLength/${BRAIN}ReducedJobs.txt

matlab -nodisplay -nosplash -r "Initalise_fetal_surface('${BRAIN}',${HEMI},${POINTSPERJOB}); exit"

fi