function Outputs = runCustomGenPhysMdl(A,A_dist,PD1,PD2,Input)

% This is a wrapper for all the physiological models in the paper.
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
%   The following are optional specifications for Input (will be set to
%   defaults if not set):
%
%   PD1Func = 'power-law' or 'exponential'. Controls the interaction between
%   PD1 and eta. Defaults to 'exponential'
%
%   PD2Func = 'power-law' or 'exponential'. Controls the interaction between
%   PD2 and lam. Defaults to 'power-law'


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


end
