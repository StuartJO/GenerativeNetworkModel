function runGenTopoMdlCV(TYPE,GROWTH,MdlNum,ITER)

% Input.Growth = GROWTH;
% Input.DistFunc = 'exponential';
% Input.IndvDist = 0;
% Input.GeneFunc = 'powerlaw';
% Input.normsum = 0;

if TYPE == 1
%    Input.AddMult = 'Mult';
    fileoutname = 'mult2';
   
elseif TYPE == 2
%    Input.AddMult = 'Add';
    fileoutname = 'add2';
    
elseif TYPE == 3
%   Input.AddMult = 'Add';
fileoutname = 'add3';
    
end

Mdlouts = load(['Revisedmdls_',fileoutname,'_Growth_',num2str(GROWTH),'_output.mat']);

Input = Mdlouts.Input{MdlNum};

addpath /projects/kg98/stuarto/BCT

mdldata = load('random200_data4topomdl.mat');

adjs = mdldata.adjs;

A_dist = mdldata.A_dist;
    
if GROWTH
    D = mdldata.dists;
else
    D = A_dist;
end   

PD = [];

P = Mdlouts.MinMdlFitParams{MdlNum};
Fcv = CrossValidateModel(adjs,A_dist,D,PD,P,1,Input);

Fcv.P = P;

output_location = '/fs02/hf49/Stuart/GrowthModel_newParc/Paper_topo_mdls/Crossvalidated/';

save([output_location,'CrossValidate_RevisedMdls_',fileoutname,'_mdl_',num2str(MdlNum),'_Growth_',num2str(GROWTH),'_iter_',num2str(ITER),'.mat'],'-struct','Fcv','-v7.3')
