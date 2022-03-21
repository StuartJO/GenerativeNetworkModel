function runGenPhysMdlCV(PARC,GROWTH,MdlNum,ITER,INPUTLOC,OUTPUTLOC)

% This is a wrapper for all the cross-validated physiological models in the
% paper.
%
% Inputs:
%
% PARC = value between 1-3 to indicate which parcellation/data to use,1 =
%   random200, 2 = Schaefer200, 3 = random200 with 18 timepoints (remove
%   adult timepoints)
%
% GROWTH = 1 to run a growth model, 0 for a static model
%
% MdlNum = a value between 1:10 indicating which model corssvalidation is
%   being performed for
%
% ITER = the iteration number (i.e., what number repeat of crossvalidation 
%   this is)
%
% INPUTLOC = location of the input data (saved output of
%   CompileGenMdlOutputs.m)
%
% OUTPUTLOC = location to save output 

if PARC == 1
    mdldata = load('random200_data4physmdl.mat');
    parcname = 'random200';
elseif PARC == 2
    mdldata = load('Schaefer200_data4physmdl.mat');
    parcname = 'Schaefer200';
elseif PARC == 3
    mdldata = load('random200_data4physmdl_18timepoints.mat');
    parcname = 'random200_18tps';
end

display(['Performing cross-validation for ',parcname,', phys growth ',num2str(GROWTH),' model number ',num2str(MdlNum),', iteration ',num2str(ITER)]);

adjs = mdldata.adjs;
A_dist = mdldata.A_dist;
  
if GROWTH
    D = mdldata.dists;
else
    D = A_dist;
end   

if MdlNum == 1

PD = [];

elseif MdlNum == 2

PD = mdldata.cCGE;

elseif MdlNum == 3

D = mdldata.cCGE;
PD = [];

elseif MdlNum == 4

PD = mdldata.uCGE;

elseif MdlNum == 5

D = mdldata.uCGE;
PD = [];

elseif MdlNum == 6

D = mdldata.hist_mpc;
PD = [];


elseif MdlNum == 7

PD = mdldata.hist_mpc;

elseif MdlNum == 8

D = mdldata.t1t2_mpc;
PD = [];

elseif MdlNum == 9

PD = mdldata.t1t2_mpc;

elseif MdlNum == 10

PD = [];

end

Mdlouts = load([INPUTLOC,'/',parcname,'_PhysMdls_Growth_',num2str(GROWTH),'_output.mat']);

Input = Mdlouts.Inputs{MdlNum};
Input.NNodes = length(A_dist);  
P = Mdlouts.OptimMdl{MdlNum}.min_maxKS.P;
Fcv = CrossValidateModel(adjs,A_dist,D,PD,P,1,Input);

Fcv.P = P;

save([OUTPUTLOC,'/CrossValidate_',parcname,'_PhysMdls_mdl_',num2str(MdlNum),'_Growth_',num2str(GROWTH),'_iter_',num2str(ITER),'.mat'],'-struct','Fcv','-v7.3')
