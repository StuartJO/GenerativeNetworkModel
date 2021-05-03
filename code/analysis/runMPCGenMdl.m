function runMPCGenMdl(SUB,GROWTH)

Input.Growth = GROWTH;
Input.ModelNum = 1;

for i = 1:4

if i == 1
    
Input.DistFunc = 'powerlaw';
Input.IndvDist = 0;
Input.GeneFunc = 'powerlaw';
Input.normsum = 0;
Input.ParamRange(1,:) = [0 50];
    
      Input.AddMult = 'Mult';
%    fileoutname = 'mult2';
% gamma
Input.ParamRange(2,:) = [0 0];
% lambda (PD2 param)
Input.ParamRange(5,:) = [0 0];
% alpha
Input.ParamRange(3,:) = [1 1];
% alpha2 (alpha for PD2)
Input.ParamRange(4,:) = [1 1];   

elseif i == 2
    Input.DistFunc = 'exponential';
Input.IndvDist = 0;
Input.GeneFunc = 'powerlaw';
Input.normsum = 0;
Input.ParamRange(1,:) = [-2 0];
    
    Input.AddMult = 'Mult';
%    fileoutname = 'mult2';
% gamma
Input.ParamRange(2,:) = [0 0];
% lambda (PD2 param)
Input.ParamRange(5,:) = [0 50];
% alpha
Input.ParamRange(3,:) = [1 1];
% alpha2 (alpha for PD2)
Input.ParamRange(4,:) = [1 1];     
elseif i == 3
    Input.DistFunc = 'exponential';
Input.IndvDist = 0;
Input.GeneFunc = 'powerlaw';
Input.normsum = 0;
Input.ParamRange(1,:) = [-2 0];
    
    Input.AddMult = 'Add';
%    fileoutname = 'add2';
% gamma
Input.ParamRange(2,:) = [1 1];
% lambda (PD2 param)
Input.ParamRange(5,:) = [0 0];
% alpha
Input.ParamRange(3,:) = [0 0];
% alpha2 (alpha for PD2)
Input.ParamRange(4,:) = [0 .05];    
    
elseif i == 4
    Input.DistFunc = 'exponential';
Input.IndvDist = 0;
Input.GeneFunc = 'powerlaw';
Input.normsum = 0;
Input.ParamRange(1,:) = [-2 0];
    
Input.AddMult = 'Add';
%fileoutname = 'add3';
% gamma
% gamma
Input.ParamRange(2,:) = [0 0];
% lambda (PD2 param)
Input.ParamRange(5,:) = [0 50];
% alpha
Input.ParamRange(3,:) = [0 0];
% alpha2 (alpha for PD2)
Input.ParamRange(4,:) = [0 8];   
    
end

addpath /projects/kg98/stuarto/BCT

%load('MPCrandom200.mat','MPCc');
load('MPCarandom200.mat','MPCa');

load('random200_100randsubs_for_modelling.mat','adjs')

A = double(adjs{SUB}>0);

load('random200_distances.mat','dists')

A_dist = dists{19};

if i == 1

    D = MPCa;

    PD = [];
    
else

    if GROWTH
        D = dists;
    else
        D = A_dist;
    end
    
    PD = MPCa;
    
end

[~,E{i},K{i},P{i},b{i},C{i},B{i},Ebest{i},Kbest{i},Cbest{i},b_best{i},EbestCorr{i},KbestCorr{i},CbestCorr{i},b_bestCorr{i}]= GenMdl(A,A_dist,D,PD,Input);

end

save(['MPCmdls_abs_Sub_',num2str(SUB),'_Growth_',num2str(GROWTH),'.mat'],'E','K','P','b','C','B','Ebest','Kbest','Cbest','b_best','EbestCorr','KbestCorr','CbestCorr','b_bestCorr')
