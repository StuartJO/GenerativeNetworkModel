#!/bin/env bash
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH -t 10:0:0
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --export=ALL
#SBATCH --mem-per-cpu=1G
#SBATCH --account=kg98
#SBATCH --array=1-100%40
#SBATCH --output=/fs02/hf49/Stuart/GrowthModel_newParc/GenerativeNetworkModel/SLURM_OUTPUT/slurm-%j_%a.out

declare -i SUB=$SLURM_ARRAY_TASK_ID
echo $SLURM_ARRAY_TASK_ID

GROWTH=$1
PARC=$2
module load matlab/r2017b

CODEDIR="/fs02/hf49/Stuart/GrowthModel_newParc/GenerativeNetworkModel"
DATADIR="/fs02/hf49/Stuart/GrowthModel_newParc/GenerativeNetworkModel/data"

matlab -nodisplay -nosplash -r "addpath(genpath(('${CODEDIR}'))); addpath('/projects/kg98/stuarto/BCT'); runGenPhysMdl($SUB,$GROWTH,$PARC,'${DIRDATA}'); exit"

