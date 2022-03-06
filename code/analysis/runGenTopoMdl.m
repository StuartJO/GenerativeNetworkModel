function runGenTopoMdl(TYPE,SUB,GROWTH,SAVELOC)

% This runs the topological generative models.
% TYPE = value between 1 to 13.
%       1.  'sptl'          spatial model
%       2.  'neighbors'     number of common neighbors
%       3.  'matching'      matching index
%       4.  'clu-avg'       average clustering coeff.
%       5.  'clu-min'       minimum clustering coeff.
%       6.  'clu-max'       maximum clustering coeff.
%       7.  'clu-diff'      difference in clustering coeff.
%       8.  'clu-prod'      product of clustering coeff.
%       9.  'deg-avg'       average degree
%       10. 'deg-min'       minimum degree
%       11. 'deg-max'       maximum degree
%       12. 'deg-diff'      difference in degree
%       13. 'deg-prod'      product of degree
%
% SUB = which subject to run the model for (1 to 100)
%
% GROWTH = set to 1 to run the growth form of the model (set to 0 otherwise) 
%
% SAVELOC = location to save the output

if nargin < 4
    SAVELOC = '.';
end


Input.useParfor = 1;

Input.ndraw = 2000;
Input.pow = 2;
Input.nlvl = 5;
Input.useParFor = 1;

Input.Growth = GROWTH;
Input.PD1Func = 'exponential';
Input.TopoFunc = 'powerlaw';
Input.seed = [];

Input.normsum = 0;

% eta (PD param)
Input.ParamRange(1,:) = [-2 0];

if TYPE == 1
    Input.AddMult = 'Mult';
    fileoutname = 'mult2';
% gamma
Input.ParamRange(2,:) = [-8 8];
% lambda (PD param)
Input.ParamRange(5,:) = [0 0];
% alpha
Input.ParamRange(3,:) = [1 1];
% alpha2 (alpha for PD)
Input.ParamRange(4,:) = [1 1];     
elseif TYPE == 2
    Input.AddMult = 'Add';
    fileoutname = 'add2';
% gamma
Input.ParamRange(2,:) = [1 1];
% lambda (PD param)
Input.ParamRange(5,:) = [0 0];
% alpha
Input.ParamRange(3,:) = [0 8];
% alpha2 (alpha for PD)
Input.ParamRange(4,:) = [0 0];    
    
elseif TYPE == 3
Input.AddMult = 'Add';
fileoutname = 'add3';
% gamma
Input.ParamRange(2,:) = [-8 8];
% lambda (PD param)
Input.ParamRange(5,:) = [0 0];
% alpha
Input.ParamRange(3,:) = [0 8];
% alpha2 (alpha for PD)
Input.ParamRange(4,:) = [0 0];   
    
end
if ~isfile([SAVELOC,'/random200_TopoMdls_Sub_',num2str(SUB),'_',fileoutname,'_Growth_',num2str(GROWTH),'.mat'])

load('random200_data4topomdl.mat','adjs','A_dist','dists')

A = double(adjs{SUB}>0);

if GROWTH
    D = dists;
else
    D = A_dist;
end

Input.NNodes = length(A_dist);

for i = 1:13
Input.ModelNum=i;
if TYPE == 2
    % Experimentation led us to find this parameter range was needed when
    % no gamma is used
    if ismember(i,[4 6 7 9 11 12 13])
        Input.ParamRange(3,:) = [0 .05];
    else
       Input.ParamRange(3,:) = [0 8]; 
    end
end
Output = GenMdl(A,A_dist,D,[],Input);

Outputs.maxKS{i} = Output.maxKS;
Outputs.DegCorr{i} = Output.DegCorr;
Outputs.KS{i} = Output.KS;
Outputs.P{i} = Output.P;
Outputs.b{i} = Output.b;

Outputs.optim_maxKS{i} = Output.optim_maxKS;
Outputs.optim_KS{i} = Output.optim_KS;
Outputs.optim_b{i} = Output.optim_b;
Outputs.optim_DegCorr{i} = Output.optim_DegCorr;

Outputs.bestDegCorr_maxKS{i} = Output.bestDegCorr_maxKS;
Outputs.bestDegCorr_KS{i} = Output.bestDegCorr_KS;
Outputs.bestDegCorr_b{i} = Output.bestDegCorr_b;
Outputs.bestDegCorr_DegCorr{i} = Output.bestDegCorr_DegCorr;

% Save the input configurations to output. Helps to keep track of what was
% done
Outputs.Input{i} = Output.Input;

end

save([SAVELOC,'/random200_TopoMdls_Sub_',num2str(SUB),'_',fileoutname,'_Growth_',num2str(GROWTH),'.mat'],'-struct','Outputs','-v7.3')

end