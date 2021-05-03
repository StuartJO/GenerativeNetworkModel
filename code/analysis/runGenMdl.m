function Outputs = runGenMdl(TYPE,SUB,GROWTH)

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

addpath /projects/kg98/stuarto/BCT

load('random200_100randsubs_for_modelling.mat','adjs')

A = double(adjs{SUB}>0);

load('random200_distances.mat','dists')

A_dist = dists{19};

if GROWTH
    D = dists;
else
    D = A_dist;
end

for i = 1:13
Input.ModelNum=i;
if TYPE == 2
    if ismember(i,[4 6 7 9 11 12 13])
        Input.ParamRange(3,:) = [0 .05];
    else
       Input.ParamRange(3,:) = [0 8]; 
    end
end
[Outputs,E{i},K{i},P{i},b{i},C{i},B{i},Ebest{i},Kbest{i},Cbest{i},b_best{i},EbestCorr{i},KbestCorr{i},CbestCorr{i},b_bestCorr{i}]= GenMdl(A,A_dist,D,[],Input);

end

save(['Revisedmdls_Sub_',num2str(SUB),'_',fileoutname,'_Growth_',num2str(GROWTH),'.mat'],'E','K','P','b','C','B','Ebest','Kbest','Cbest','b_best','EbestCorr','KbestCorr','CbestCorr','b_bestCorr')
