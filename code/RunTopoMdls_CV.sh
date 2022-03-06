#!/bin/env bash
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH -t 3-0:0:0
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --export=ALL
#SBATCH --mem-per-cpu=10G
#SBATCH --account=kg98
#SBATCH --output=/fs02/hf49/Stuart/GrowthModel_newParc/GenerativeNetworkModel/SLURM_OUTPUT/slurm-%j_%a.out
#SBATCH --array=1-20%20

declare -i ITER=$SLURM_ARRAY_TASK_ID
echo $SLURM_ARRAY_TASK_ID

#for mdl in 3 1 2 4 5 6 7 8 9 10 11 12 13; do for type in 1 2 3; do for growth in 0 1; do sbatch ./RunTopoMdls_CV.sh $type $growth $mdl; done; done; done

#SUB=$1
TYPE=$1
GROWTH=$2
MDLNUM=$3

echo $1 $2
module load matlab/r2017a
#module load matlab

matlab -nodisplay -nosplash -r "addpath(genpath(('/fs02/hf49/Stuart/GrowthModel_newParc/GenerativeNetworkModel'))); addpath('/projects/kg98/stuarto/BCT'); runGenTopoMdlCV($TYPE,$GROWTH,$MDLNUM,$ITER,'/fs02/hf49/Stuart/GrowthModel_newParc/GenerativeNetworkModel/data/OUTPUTS','/fs02/hf49/Stuart/GrowthModel_newParc/GenerativeNetworkModel/data/Crossvalidated/MULT2'); exit"
