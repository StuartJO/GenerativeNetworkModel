function runGenTopoMdlCV(TYPE,GROWTH,MdlNum,ITER,INPUTLOC,OUTPUTLOC)

%OUTPUTLOC = '/fs02/hf49/Stuart/GrowthModel_newParc/Paper_topo_mdls/Crossvalidated/';

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

if ~isfile([OUTPUTLOC,'/CrossValidate_random200_TopoMdls_',fileoutname,'_mdl_',num2str(MdlNum),'_Growth_',num2str(GROWTH),'_iter_',num2str(ITER),'.mat'])

Mdlouts = load([INPUTLOC,'/random200_TopoMdls_',fileoutname,'_Growth_',num2str(GROWTH),'_output.mat'],'Inputs','OptimMdl');

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

Input = Mdlouts.Inputs{MdlNum};
Input.NNodes = length(A_dist);  
P = Mdlouts.OptimMdl{MdlNum}.min_maxKS.P;
Fcv = CrossValidateModel(adjs,A_dist,D,PD,P,1,Input);

Fcv.P = P;

save([OUTPUTLOC,'/CrossValidate_random200_TopoMdls_',fileoutname,'_mdl_',num2str(MdlNum),'_Growth_',num2str(GROWTH),'_iter_',num2str(ITER),'.mat'],'-struct','Fcv','-v7.3')

end