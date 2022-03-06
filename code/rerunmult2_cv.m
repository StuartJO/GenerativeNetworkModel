addpath(genpath(('/fs02/hf49/Stuart/GrowthModel_newParc/GenerativeNetworkModel'))); 
addpath('/projects/kg98/stuarto/BCT'); 

for i = 13:-1:4
    parfor j = 1:20
        runGenTopoMdlCV(1,0,i,j,'/fs02/hf49/Stuart/GrowthModel_newParc/GenerativeNetworkModel/data/OUTPUTS','/fs02/hf49/Stuart/GrowthModel_newParc/GenerativeNetworkModel/data/Crossvalidated/MULT2')
    end
end