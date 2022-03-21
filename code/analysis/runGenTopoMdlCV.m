function runGenTopoMdlCV(TYPE,GROWTH,MdlNum,ITER,INPUTLOC,OUTPUTLOC)

% This is a wrapper for all the cross-validated topological models in the
% paper.
%
% Inputs:
%
% TYPE = value between 1 to 3.
%       1.  'mult2'    multiplicative model
%       2.  'add2'     additive model without a gamma parameter
%       3.  'add3'     additive model with a gamma parameter
%
% GROWTH = 1 to run a growth model, 0 for a static model
%
% MdlNum = a value between 1:10 indicating which model corssvalidation is
%   being performed for
%
% ITER = the iteration number (i.e., what number repeat of crossvalidation 
%   this is)
%
%
% INPUTLOC = location of the input data (saved output of
%   CompileGenMdlOutputs.m)
%
% OUTPUTLOC = location to save output 

if TYPE == 1

    fileoutname = 'mult2';
   
elseif TYPE == 2

    fileoutname = 'add2';
    
elseif TYPE == 3

    fileoutname = 'add3';
    
end

if ~isfile([OUTPUTLOC,'/CrossValidate_random200_TopoMdls_',fileoutname,'_mdl_',num2str(MdlNum),'_Growth_',num2str(GROWTH),'_iter_',num2str(ITER),'.mat'])

Mdlouts = load([INPUTLOC,'/random200_TopoMdls_',fileoutname,'_Growth_',num2str(GROWTH),'_output.mat'],'Inputs','OptimMdl');

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

else
   disp([OUTPUTLOC,'/CrossValidate_random200_TopoMdls_',fileoutname,'_mdl_',num2str(MdlNum),'_Growth_',num2str(GROWTH),'_iter_',num2str(ITER),'.mat exists!'] )
end