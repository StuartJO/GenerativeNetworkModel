function runGenTopoMdl(TYPE,SUB,GROWTH,SAVELOC,MDLS,save_suffix)

% This runs the topological generative models.
% TYPE = value between 1 to 3.
%       1.  'mult2'    multiplicative model
%       2.  'add2'     additive model without a gamma parameter
%       3.  'add3'     additive model with a gamma parameter
%
% SUB = which subject to run the model for (1 to 100)
%
% GROWTH = set to 1 to run the growth form of the model (set to 0 otherwise) 
%
% SAVELOC = location to save the output
%
% MDLS = The index of topological models to run

if nargin < 4
    SAVELOC = '.';
end

if nargin < 5
    MDLS = 1:13;
end

if nargin < 6
    save_suffix = '';
end

Input.useParfor = 1;

Input.ndraw = 2000;
Input.pow = 2;
Input.nlvl = 5;

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
Input.ParamRange(5,:) = [NaN NaN];
% alpha
Input.ParamRange(3,:) = [NaN NaN];
% alpha2 (alpha for PD)
Input.ParamRange(4,:) = [NaN NaN];     
elseif TYPE == 2
    Input.AddMult = 'Add';
    fileoutname = 'add2';
% gamma
Input.ParamRange(2,:) = [1 1];
% lambda (PD param)
Input.ParamRange(5,:) = [NaN NaN];
% alpha
Input.ParamRange(3,:) = [0 8];
% alpha2 (alpha for PD)
Input.ParamRange(4,:) = [NaN NaN];    
    
elseif TYPE == 3
Input.AddMult = 'Add';
fileoutname = 'add3';
% gamma
Input.ParamRange(2,:) = [-8 8];
% lambda (PD param)
Input.ParamRange(5,:) = [NaN NaN];
% alpha
Input.ParamRange(3,:) = [0 8];
% alpha2 (alpha for PD)
Input.ParamRange(4,:) = [NaN NaN];   
    
end

load('random200_data4topomdl.mat','adjs','A_dist','dists')

A = double(adjs{SUB}>0);

if GROWTH
    D = dists;
else
    D = A_dist;
end

Input.NNodes = length(A_dist);

MDLIND = 1;

for MDL = MDLS
    
Input.ModelNum=MDL;

if MDL == 1
    Input.ParamRange(2:5,:) = NaN;
end

if TYPE == 2
    % Experimentation led us to find this parameter range was needed when
    % no gamma is used
    if ismember(MDL,[4 6 7 9 11 12 13])
        Input.ParamRange(3,:) = [0 .05];
    else
       Input.ParamRange(3,:) = [0 8]; 
    end
end
Output = GenMdl(A,A_dist,D,[],Input);

Outputs.maxKS{MDLIND} = Output.maxKS;
Outputs.DegCorr{MDLIND} = Output.DegCorr;
Outputs.KS{MDLIND} = Output.KS;
Outputs.P{MDLIND} = Output.P;
Outputs.b{MDLIND} = Output.b;

Outputs.optim_maxKS{MDLIND} = Output.optim_maxKS;
Outputs.optim_KS{MDLIND} = Output.optim_KS;
Outputs.optim_b{MDLIND} = Output.optim_b;
Outputs.optim_DegCorr{MDLIND} = Output.optim_DegCorr;

Outputs.bestDegCorr_maxKS{MDLIND} = Output.bestDegCorr_maxKS;
Outputs.bestDegCorr_KS{MDLIND} = Output.bestDegCorr_KS;
Outputs.bestDegCorr_b{MDLIND} = Output.bestDegCorr_b;
Outputs.bestDegCorr_DegCorr{MDLIND} = Output.bestDegCorr_DegCorr;

% Save the input configurations to output. Helps to keep track of what was
% done
Outputs.Input{MDLIND} = Output.Input;

MDLIND = MDLIND + 1;

end

save([SAVELOC,'/random200_TopoMdls_Sub_',num2str(SUB),'_',fileoutname,'_Growth_',num2str(GROWTH),save_suffix,'.mat'],'-struct','Outputs','-v7.3')
