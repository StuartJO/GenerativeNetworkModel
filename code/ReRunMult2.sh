#!/bin/env bash

module load matlab

Nfiles=$(ls -1U /fs02/hf49/Stuart/GrowthModel_newParc/GenerativeNetworkModel/data/Optimisation/MULT2/ | wc -l)

while [ $Nfiles -lt 100 ]; do
	sleep 60
	Nfiles=$(ls -1U /fs02/hf49/Stuart/GrowthModel_newParc/GenerativeNetworkModel/data/Optimisation/MULT2/ | wc -l)
	echo $Nfiles
done
#sleep 5m
echo "Files generated!"
matlab -nodisplay -nosplash -r "addpath(genpath(('/fs02/hf49/Stuart/GrowthModel_newParc/GenerativeNetworkModel'))); addpath('/projects/kg98/stuarto/BCT'); GetNEWTopoMdlMultResult; exit"

for mdl in 3 1 2 4 5 6 7 8 9 10 11 12 13; do 
	for type in 1; do 
		for growth in 0; do 
			sbatch ./RunTopoMdls_CV.sh $type $growth $mdl 
		done 
	done 
done

