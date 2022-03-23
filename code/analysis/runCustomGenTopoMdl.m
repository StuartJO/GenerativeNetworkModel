function MdlOutput = runCustomGenTopoMdl(A,A_dist,D,TYPE,MDLS)

% This runs the topological generative models with the form
% (PD1^eta)*a1(T^gam)or (PD1^eta)+a1(T^gam)
%
% where PD is a measure of distance/similarity between a pair of
% nodes; T is some measure of the topology between a pair of nodes; and
% eta, gam, and a1 are free parameters. Note that each power-law
% interaction (e.g., PD1^eta) can be replaced with an exponential function
% (e.g., exp(eta*PD1)).
%
% Inputs:
% A = The adjacency matrix to model
%
% A_dist = The distances between nodes (e.g., Euclidean distances, fibre
% distance) for A. This distance is used to compute how similar the model
% is to the input adjacency matrix (A).
%
% D = Either a matrix indicating pairwise similarity/distances between nodes, or a cell where each
% element contains such a distance matrix, to be used for the modelling itself.
% The "growth" model is specified by the latter of these options.
%
% TYPE = a string with the following format
%       'mult2'    multiplicative model
%       'add2'     additive model without a gamma parameter
%       'add3'     additive model with a gamma parameter
% MDLS = The index of topological models to run corresponding to the
% following:
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
% Outputs:
% MdlOutput = a structure with the following fields:
%   maxKS{MDL} = N*1 array of maxKS values for each network generated
%       during optimisation for model number 'MDL'
%   DegCorr{MDL} = N*1 array of values for the correlation with  
%       empirical degree each network generated during optimisation for model number 'MDL'  
%   KS{MDL} = N*4 matrix of KS values for each network generated
%       during optimisation for model number 'MDL'. Each column corresponds to the following 
%       measure: 1 = degree; 2 = clustering; 3 = betweenness; 4 = mean edge 
%       length
%   P{MDL} = N*5 matrix of parameter values for each network generated
%       during optimisation for model number 'MDL'
%   b{MDL} = 1*N cell array, where each cell contains an edge index
%       list for each network generated during optimisation for model number 'MDL'
% 
%   optim_maxKS{MDL} = for networks generated using the parameters
%       of the lowest maxKS value for model number 'MDL', a 1*100 array of maxKS values
%   optim_KS{MDL} = for networks generated using the parameters
%       of the lowest maxKS value for model number 'MDL', an 100*4 matrix of KS values
%   optim_b{MDL} = a 1*100 cell array where each cell contains an edge
%       index list for each network generated during using the parameters 
%       of the lowest maxKS value for model number 'MDL'       
%   optim_DegCorr{MDL} = for networks generated using the parameters
%       of the lowest maxKS value for model number 'MDL', a 1*100 array of values for the 
%       correlation with empirical degree      
%
%   bestDegCorr_maxKS{MDL} = for networks generated using the parameters
%       of the largest degree correlation for model number 'MDL', a 1*100 array of maxKS values
%   bestDegCorr_KS{MDL} = for networks generated using the parameters
%       of the largest degree correlation for model number 'MDL', an 100*4 matrix of KS values
%   bestDegCorr_b{MDL} = a 1*100 cell array where each cell contains an edge
%       index list for each network generated during using the parameters 
%       of the largest degree correlation for model number 'MDL'    
%   bestDegCorr_DegCorr{MDL} = for networks generated using the parameters
%       of the largest degree correlation for model number 'MDL', a 1*100 array of values for the 
%       correlation with empirical degree     
if nargin < 5
    MDLS = 1:13;
end

Input.useParfor = 1;

Input.ndraw = 2000;
Input.pow = 2;
Input.nlvl = 5;

Input.PD1Func = 'exponential';
Input.TopoFunc = 'powerlaw';
Input.seed = [];

Input.normsum = 0;

% eta (PD param)
Input.ParamRange(1,:) = [-2 0];

if strcmp(TYPE,'mult2')
    Input.AddMult = 'Mult';
% gamma
Input.ParamRange(2,:) = [-8 8];
% lambda (PD param)
Input.ParamRange(5,:) = [NaN NaN];
% alpha
Input.ParamRange(3,:) = [NaN NaN];
% alpha2 (alpha for PD)
Input.ParamRange(4,:) = [NaN NaN];     
elseif strcmp(TYPE,'add2')
    Input.AddMult = 'Add';
% gamma
Input.ParamRange(2,:) = [NaN NaN];
% lambda (PD param)
Input.ParamRange(5,:) = [NaN NaN];
% alpha
Input.ParamRange(3,:) = [0 8];
% alpha2 (alpha for PD)
Input.ParamRange(4,:) = [NaN NaN];    
    
elseif strcmp(TYPE,'add3')
Input.AddMult = 'Add';
% gamma
Input.ParamRange(2,:) = [-8 8];
% lambda (PD param)
Input.ParamRange(5,:) = [NaN NaN];
% alpha
Input.ParamRange(3,:) = [0 8];
% alpha2 (alpha for PD)
Input.ParamRange(4,:) = [NaN NaN];   
else
    error('Unknown ''TYPE'' specified')
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

MdlOutput.maxKS{MDLIND} = Output.maxKS;
MdlOutput.DegCorr{MDLIND} = Output.DegCorr;
MdlOutput.KS{MDLIND} = Output.KS;
MdlOutput.P{MDLIND} = Output.P;
MdlOutput.b{MDLIND} = Output.b;

MdlOutput.optim_maxKS{MDLIND} = Output.optim_maxKS;
MdlOutput.optim_KS{MDLIND} = Output.optim_KS;
MdlOutput.optim_b{MDLIND} = Output.optim_b;
MdlOutput.optim_DegCorr{MDLIND} = Output.optim_DegCorr;

MdlOutput.bestDegCorr_maxKS{MDLIND} = Output.bestDegCorr_maxKS;
MdlOutput.bestDegCorr_KS{MDLIND} = Output.bestDegCorr_KS;
MdlOutput.bestDegCorr_b{MDLIND} = Output.bestDegCorr_b;
MdlOutput.bestDegCorr_DegCorr{MDLIND} = Output.bestDegCorr_DegCorr;

% Save the input configurations to output. Helps to keep track of what was
% done
MdlOutput.Input{MDLIND} = Output.Input;

MDLIND = MDLIND + 1;

end