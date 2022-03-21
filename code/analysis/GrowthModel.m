function [B,b] = GrowthModel(PD1,PD2,P,m,Input)

% This function will generate a network while accounting for changes in the
% distances of nodes over time.
% The model has the following basic form for defining connection probabilities:
%
% (PD1^eta)*a1(T^gam)*a2(PD2^lam) or (PD1^eta)+a1(T^gam)+a2(PD^lam)
%
% where PD1/PD2 is a measure of distance/similarity between a pair of
% nodes; T is some measure of the topology between a pair of nodes; and
% eta, gam, lam, a1, and a2 are free parameters. Note that each power-law
% interaction (e.g., PD1^eta) can be replaced with an exponential function
% (e.g., exp(eta*PD1)).
%
% See the bottom of this header as to the model form appears different to 
% what is in the paper/why there are seemingly more than 3 free parameters
%
% Teachnically, this function is a wrapper for another function
% (GrowthModel) which itself is a wrapper for functions which actually run
% the generative model. Turtles all the way down.
%
% Inputs:
%
% PD1 = Either a matrix indicating pairwise similarity/distances between nodes, or a cell where each
% element contains such a distance matrix, to be used for the modelling itself.
% The "growth" model is specied by the latter of these options.
%
% PD2 = A second matrix indicating pairwise similarity/distances between nodes in
% A (unlike with PD1 a cell cannot be used as an input here). This could be 
% correlated gene expression, similarity in histology etc.
%
% m = the number of edges for the model to form. When doing a growth model
% this should be a vector where each elemenet specifies the number of edges
% which should be present in the network at the end of that timestep (it is
% NOT the number of edges to form at each step, but rather the cumulative
% edges). Can also be specified as the density value
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
%   PD1Func = 'power-law' or 'exponential'. Controls the interaction between
%   PD1 and eta. Defaults to 'exponential'
%
%   TopoFunc = 'power-law' or 'exponential'. Controls the interaction between T
%   and eta. Defaults to 'power-law'
%
%   The following are optional specifications for Input (will be set to
%   defaults if not set):
%
%       PD2Func = 'power-law' or 'exponential'. Controls the interaction between
%       PD2 and lam. Defaults to 'power-law'. If PD2 ~= [] then this is
%       required
%
%       epsilon = an amount to add to each edges topology value in the model (to
%       ensure each edge doesn't become undefinied). Defaults to 0 when 
%       Input.AddMult = 'Add' or 1e-6 when Input.AddMult = 'Mult'
%
%       seed = a seed network. If none is desired set to [] (default)
%
%       normsum = set to 1 to normalise each term by its sum (set to 0 by
%       default, which normalises by the max)
%
% Outputs:
%
% B = an adjanceny matrix of the generated network
%
% b = an array of edge indices 
%
% Written by Stuart Oldham, 2021

%% Check inputs
if iscell(PD1)
    N = length(PD1{1});
    Dists = PD1;
else 
    N = length(PD1);
    Dists{1} = PD1;
end

if iscell(PD1)
    if length(PD1) ~= length(m)  
        error('PD1 and m must be of the same length')
    end
end

if ~isfield(Input,'epsilon')
    switch Input.AddMult
    case 'Add'
    Input.epsilon = 0;
    case 'Mult'
    Input.epsilon = 1e-6;
    end
end

if ~isfield(Input,'normsum')
    normsum = 0;  
else
    normsum = Input.normsum;  
end

TopoType = Input.TopoType;

epsilon = Input.epsilon;

AddMult = Input.AddMult;

if isempty(Input.seed) || ~isfield(Input,'seed')
   Seed = zeros(N);
else
   Seed = Input.seed; 
end

if ~isfield(Input,'PD2Func')
    if ~isempty(PD2)
        error('PD2 exists but PD2Func is not defined')
    end
    % If not defined we just set it to 'powerlaw', won't do anything
    PD2Func = 'powerlaw';
else
    PD2Func = Input.PD2Func;
end

if ~iscell(PD2) && ~isempty(PD2)
    PD2 = num2cell(PD2,[1 2]);
end

modelvar = cell(1,3);
modelvar{3} = PD2Func;
%% Set up the functions for the generative rule

PD1Func = Input.PD1Func;
TopoFunc = Input.TopoFunc;

if strcmp(PD1Func,'exponential')
    modelvar{1} = 'exponential';
elseif strcmp(PD1Func,'powerlaw')
    modelvar{1} = 'powerlaw';
else
    error('Unknown input for PD1Func')
end

if strcmp(TopoFunc,'exponential')
    modelvar{2} = 'exponential';
elseif strcmp(TopoFunc,'powerlaw')
    modelvar{2} = 'powerlaw';
else
    error('Unknown input for TopoFunc')
end

if ~isempty(PD2)
    if strcmp(PD2Func,'exponential')
        modelvar{3} = 'exponential';
    elseif strcmp(PD2Func,'powerlaw')
        modelvar{3} = 'powerlaw';
    else
        error('Unknown input for PD2Func')
    end
end

%% Run generative model

% Create network

B = Seed;

eta = P(1);
gam = P(2);
a1 = P(3);
a2 = P(4);
lam = P(5);

PDexpo = [eta lam];

% Three alpha_vals are used, however we only talk about 2 in the paper.
% Simply we set up the code to have an alpha value attached to PD1. However
% we subsequently realised that this was redundant because due to the
% normalisation, only a single alpha value is needed for the other term
% (either T or PD2) to achieve the same trade-off as if two were used.
% However I never updated the code to reflect this so one is included.

alpha_vals = [1 a1 a2];

% Basically for the growth model, we just loop over the distances in the
% cell. For each iteration, we input the output matrix from the last
% iteration as the seed for the model. Repeat until done.

for I = 1:length(Dists)
        
    if m(I) <= 1
        desiredEdges = ceil(m(I)*N*(N-1)/2);
    else
        desiredEdges = m(I);
    end 
            % PD1 and PD2 are combined into one cell
            PD = [Dists{I} PD2];

            if ~iscell(PD) && ~isempty(PD)
                PD = num2cell(PD,[1 2]);
            end

            switch AddMult
                case 'Add'

                    if normsum == 0
            [B,btemp] = gen_model_add_normmax(B,PD,desiredEdges,TopoType,modelvar,PDexpo,gam,alpha_vals,epsilon);
                    else
            [B,btemp] = gen_model_add_normsum(B,PD,desiredEdges,TopoType,modelvar,PDexpo,gam,alpha_vals,epsilon);

                    end
                case 'Mult'
            [B,btemp] = gen_model_mult(B,PD,desiredEdges,TopoType,{{modelvar{1},modelvar{3}},modelvar{2}},PDexpo,gam,epsilon);             

            end

    if I == 1 
        b = btemp;
    else
            if m(I-1) <= 1
                    PreviousDesiredEdges = ceil(m(I-1)*N*(N-1)/2);
            else
                    PreviousDesiredEdges = m(I-1);
            end 

        b = [b; btemp(PreviousDesiredEdges+1:desiredEdges)];
    end
                    

end

end