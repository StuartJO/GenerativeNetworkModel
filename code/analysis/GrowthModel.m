function [B,b] = GrowthModel(PD1,PD2,P,m,Input)

%(AddOrMult,D,TopoType,Eta,Gamma,alpha,Kalpha,m,AddSingle,epsilon,InitialSEED,DistFunc,TopoFunc,PDM,PDMFunc,normsum)

% This function will generate a network while accounting for changes in the
% distances of nodes over time.

% Inputs:                           D = a distance matrix or a cell of 
%                                       distances matrices
%                            TopoType = the generative rule (see below)
%                                 Eta = the parameter controlling distance
%                               Gamma = the parameter controlling topology
%                               alpha = a vector of alpha values. alpha(1)
%                                       is used for the distance function,
%                                       alpha(2) is for the topology, while
%                                       alpha(3) is for the "O" function
%                                       (if "O" varies across X timepoints,
%                                       alpha(3:X+3) needs to be specified
%                         m = a scaler indicating the desired
%                                       density/edges of the network, or a 
%                                       vector of values indicating the 
%                                       density/edges at each desired
%                                       timepoint (must be the same length
%                                       as D)
%                           AddSingle = if set to 0, when generating a
%                                       network, at each timepoint the 
%                                       connection probabilities will be
%                                       calculated once and connections
%                                       will be added to the network
%                                       simultaneously based on these
%                                       probabilities. If set to 1,
%                                       connections are added in one at a
%                                       time and the probabilities are also
%                                       updated each time
%                                   E = the baseline probability of forming 
%                                       a connection (default = 1e-5)
%                         InitialSEED = an initial network a seed
%                                       connections. Set as empty if there
%                                       is no seed network (default)
%                            DistFunc = specifies whether the generative 
%                                       rules for distance are based on
%                                       power-law ('powerlaw') or 
%                                       exponential ('exponential') 
%                                       function (default = 'exponential')
%                            TopoFunc = specifies whether the generative 
%                                       rules for topology are based on
%                                       power-law ('powerlaw') or 
%                                       exponential ('exponential') 
%                                       function (default = 'powerlaw')
% 
% List of generative rules impliamented (adapted from code made available
% by Rick Betzel)
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
% Written by Stuart Oldham, 2018

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

TopoType = Input.TopoType;

epsilon = Input.epsilon;

AddMult = Input.AddMult;

if isempty(Input.seed) 
    Seed = zeros(N);
end

if ~isfield(Input,'PD2Func')
    if ~isempty(PD2)
        error('PD2 exists but PD2Func is not defined')
    end
    % If not defined we just set it to 'powerlaw', won't do anything
PD2Func = '';  
else
PD2Func = Input.PD2Func;
end

normsum = Input.normsum;

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