function MdlOutput = runCustomGenPhysMdl(A,A_dist,PD1,PD2,Input)

% This runs physiological models with the form
% (PD1^eta)*a1(T^gam)*a2(PD2^lam) or (PD1^eta)+a1(T^gam)+a2(PD^lam)
%
% where PD1/PD2 is a measure of distance/similarity between a pair of
% nodes; T is some measure of the topology between a pair of nodes; and
% eta, gam, lam, a1, and a2 are free parameters. Note that each power-law
% interaction (e.g., PD1^eta) can be replaced with an exponential function
% (e.g., exp(eta*PD1)).
%
% Inputs:
%
% A = The adjacency matrix to model
%
% A_dist = The distances between nodes (e.g., Euclidean distances, fibre
% distance) for A. This distance is used to compute how similar the model
% is to the input adjacency matrix (A).
%
% PD1 = Either a matrix indicating pairwise similarity/distances between nodes, or a cell where each
% element contains such a distance matrix, to be used for the modelling itself.
% The "growth" model is specified by the latter of these options.
%
% PD2 = A second matrix indicating pairwise similarity/distances between nodes in
% A (unlike with PD1 a cell cannot be used as an input here). This could be 
% correlated gene expression, similarity in histology etc.
%
% Input = a structure containing the many possible options needed to
% configure the model. The following fields are required:
%   AddMult = 'Add' or 'Mult', if the multiplicative or additive form is to
%   be used
%
%   ModelNum = A value between 1 and 13, each corresponds to the following
%   topology form:
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
%       If not specified, Input.TopoType will be the name of the topological form 
%
%  ParamRange = A 5x2 matrix where the first column gives the lower, and
%  the second the upper, bounds for the following parameters (see the basic
%  form outlined at the top):
%       ParamRange(1,:) = eta
%       ParamRange(2,:) = gam
%       ParamRange(3,:) = a1
%       ParamRange(4,:) = a2
%       ParamRange(5,:) = lam
%
%  'ParamRange(i,:) = [x x]' will mean 'x' is always used as an input for
%  parameter 'i'. 'ParamRange(i,:) = [NaN NaN]' means that parameter isn't
%  being included (if this parameter is required to be a value, you will
%  encounter an error)
%
%   The following are optional specifications for Input (will be set to
%   defaults if not set):
%
%   PD1Func = 'power-law' or 'exponential'. Controls the interaction between
%   PD1 and eta. Defaults to 'exponential'
%
%   PD2Func = 'power-law' or 'exponential'. Controls the interaction between
%   PD2 and lam. Defaults to 'power-law'
%
% Output:
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


Input.useParfor = 1;

Input.TopoFunc = 'powerlaw';
Input.seed = [];

Input.normsum = 0;
Input.ndraw = 2000;
Input.pow = 2;
Input.nlvl = 5;

MDLIND = 1;
Output = GenMdl(A,A_dist,PD1,PD2,Input);

% Things are being saved as a cell even though there is only a single
% instance because CompiileGenMdlOutputs.m assumes this is the case
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


end
