#!/bin/env bash
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH -t 3-0:0:0
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

TYPE=$1
GROWTH=$2
module load matlab/r2017b
#module load matlab

#matlab -nodisplay -nosplash -r "RunCGETest($SUB); exit"
#matlab -nodisplay -nosplash -r "RunMdlsGeneData($SUB); exit"
matlab -nodisplay -nosplash -r "addpath(genpath(('/fs02/hf49/Stuart/GrowthModel_newParc/GenerativeNetworkModel'))); addpath('/projects/kg98/stuarto/BCT'); runGenTopoMdl($TYPE,$SUB,$GROWTH,'/fs02/hf49/Stuart/GrowthModel_newParc/GenerativeNetworkModel/data/Optimisation/MULT2'); exit"

