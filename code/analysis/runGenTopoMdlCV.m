function runGenMdlCV(TYPE,GROWTH,MdlNum,ITER)

Input.Growth = GROWTH;
Input.DistFunc = 'exponential';
Input.IndvDist = 0;
Input.GeneFunc = 'powerlaw';

Input.normsum = 0;

% eta (PD1 param)
Input.ParamRange(1,:) = [-2 0];

if TYPE == 1
    Input.AddMult = 'Mult';
    fileoutname = 'mult2';
% gamma
Input.ParamRange(2,:) = [-8 8];
% lambda (PD2 param)
Input.ParamRange(5,:) = [0 0];
% alpha
Input.ParamRange(3,:) = [1 1];
% alpha2 (alpha for PD2)
Input.ParamRange(4,:) = [1 1];     
elseif TYPE == 2
    Input.AddMult = 'Add';
    fileoutname = 'add2';
% gamma
Input.ParamRange(2,:) = [1 1];
% lambda (PD2 param)
Input.ParamRange(5,:) = [0 0];
% alpha
Input.ParamRange(3,:) = [0 8];
% alpha2 (alpha for PD2)
Input.ParamRange(4,:) = [0 0];    
    
elseif TYPE == 3
Input.AddMult = 'Add';
fileoutname = 'add3';
% gamma
Input.ParamRange(2,:) = [-8 8];
% lambda (PD2 param)
Input.ParamRange(5,:) = [0 0];
% alpha
Input.ParamRange(3,:) = [0 8];
% alpha2 (alpha for PD2)
Input.ParamRange(4,:) = [0 0];   
    
end

%mdlpath = '/fs02/hf49/Stuart/GrowthModel_newParc/Paper_topo_mdls/Optimisation/';

%Mdlouts = CompileGenMdl([mdlpath,'Revisedmdls_Sub_$_',fileoutname,'_Growth_',num2str(GROWTH),'.mat'],100,[],0);

Mdlouts = load(['Revisedmdls_',fileoutname,'_Growth_',num2str(GROWTH),'_output.mat']);

addpath /projects/kg98/stuarto/BCT

load('random200_100randsubs_for_modelling.mat','adjs')

%A = double(adjs{SUB}>0);

load('random200_distances.mat','ADJS')

A_dist = ADJS{19};

if GROWTH
    D = ADJS;
else
    D = A_dist;
end

PD = [];

Input.ModelNum=MdlNum;
if TYPE == 2
    if ismember(MdlNum,[4 6 7 9 11 12 13])
        Input.ParamRange(3,:) = [0 .05];
    else
       Input.ParamRange(3,:) = [0 8]; 
    end
end

P = Mdlouts.MinMdlFitParams{MdlNum};
Fcv = CrossValidateModel(adjs,A_dist,D,PD,P,1,Input);

Fcv.P = P;

output_location = '/fs02/hf49/Stuart/GrowthModel_newParc/Paper_topo_mdls/Crossvalidated/';

save([output_location,'CrossValidate_RevisedMdls_',fileoutname,'_mdl_',num2str(MdlNum),'_Growth_',num2str(GROWTH),'_iter_',num2str(ITER),'.mat'],'-struct','Fcv','-v7.3')
