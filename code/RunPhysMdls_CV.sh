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

PARC=$1
GROWTH=$2
MDLNUM=$3

CODEDIR="/fs02/hf49/Stuart/GrowthModel_newParc/GenerativeNetworkModel"
DATADIR="/fs02/hf49/Stuart/GrowthModel_newParc/GenerativeNetworkModel/data"

echo $1 $2 $3
module load matlab/r2017a

matlab -nodisplay -nosplash -r "addpath(genpath(('${CODEDIR}'))); addpath('/projects/kg98/stuarto/BCT'); runGenPhysMdlCV($PARC,$GROWTH,$MDLNUM,$ITER,'${DATADIR}/OUTPUTS','${DATADIR}/Crossvalidated'); exit"
