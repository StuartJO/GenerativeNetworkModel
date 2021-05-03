function [B,b] = GrowthModel(AddOrMult,D,TopoType,Eta,Gamma,alpha,Kalpha,density_val,AddSingle,E,InitialSEED,DistFunc,TopoFunc,PDM,PDMFunc,normsum)

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
%                         density_val = a scaler indicating the desired
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
%                                 PDM = A pairwise distance matrix, the
%                                       same size as the matrices in D
%
%                             PDMFunc = specifies whether the generative 
%                                       rules for the PDM are based on a
%                                       power-law ('powerlaw') or 
%                                       exponential ('exponential') 
%                                       function (default = 'powerlaw') 
%
%                             normsum = set to 1 to normalise the funciton
%                                       by the sum of values rather than
%                                       the max
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
%       14. 'com'           network communicability
%
% Written by Stuart Oldham, 2018

%% Check inputs
if iscell(D)
    N = length(D{1});
    Dists = D;
else 
    N = length(D);
    Dists{1} = D;
end

if iscell(D)
    if length(D) ~= length(density_val)  
        error('D and density_val must be of the same length')
    end
end

if nargin < 9
   AddSingle = 0; 
end

if nargin < 10
    E = 1e-5;
end

if nargin < 11
    InitialSEED = zeros(N);
end

if nargin < 12
    DistFunc = 'powerlaw';
end

if nargin < 13
    TopoFunc = 'powerlaw';
end

if nargin < 14
    PDM = [];
end

if nargin < 15
    PDMFunc = 'powerlaw';
end

if nargin < 16
    normsum = 0;
end

if isempty(InitialSEED)
    InitialSEED = zeros(N);
end

if ~iscell(PDM) && ~isempty(PDM)
    PDM = num2cell(PDM,[1 2]);
end

modelvar = cell(1,3);
modelvar{3} = PDMFunc;
%% Set up the functions for the generative rule

if strcmp(DistFunc,'exponential')
    fun1 = @(x,y) exp(x.*y);
    modelvar{1} = 'exponential';
elseif strcmp(DistFunc,'powerlaw')
    fun1 = @(x,y) x.^y;
    modelvar{1} = 'powerlaw';
else
    error('Unknown input for DistFunc')
end

if strcmp(TopoFunc,'exponential')
    fun2 = @(x) exp((TopologicalFunction(x,TopoType,E).*Gamma));
    modelvar{2} = 'exponential';
elseif strcmp(TopoFunc,'powerlaw')
    fun2 = @(x) (TopologicalFunction(x,TopoType,E).^Gamma);
    modelvar{2} = 'powerlaw';
else
    error('Unknown input for TopoFunc')
end

%% Run generative model

% Create network

B = InitialSEED;

for I = 1:length(Dists)
        
        if density_val(I) <= 1
            desiredEdges = ceil(density_val(I)*N*(N-1)/2);
        else
            desiredEdges = density_val(I);
        end 
        
        if AddSingle

                PD = [Dists{I} PDM];
                
                if ~iscell(PD) && ~isempty(PD)
                    PD = num2cell(PD,[1 2]);
                end
                
                PDexpo = Eta;
                
                switch AddOrMult
                    case 'Add'
            
                        
                        if normsum == 0
                [B,btemp] = gen_model_additive_fast2_PD(B,PD,desiredEdges,TopoType,modelvar,PDexpo,Gamma,[1 Kalpha alpha(2)],E);
                        else
                 [B,btemp] = gen_model_additive_fast2_PD_normsum(B,PD,desiredEdges,TopoType,modelvar,PDexpo,Gamma,[1 Kalpha alpha(2)],E);
                           
                        end

                    case 'Mult'
                [B,btemp] = gen_model_mult_pd(B,PD,desiredEdges,TopoType,{{modelvar{1},modelvar{3}},modelvar{2}},PDexpo,Gamma,E);             
         
                end
                
        if I == 1 
			b = btemp;
		else
        		if density_val(I-1) <= 1
            			PreviousDesiredEdges = ceil(density_val(I-1)*N*(N-1)/2);
        		else
            			PreviousDesiredEdges = density_val(I-1);
        		end 

			b = [b; btemp(PreviousDesiredEdges+1:desiredEdges)];
		end
                
                
        else
            
            if strcmp('GenType','sptl')
               SEED = alpha(1)*fun1(Dists{I}); 
            else
               SEED = alpha(1)*fun1(Dists{I}) + Kalpha*fun2(B);
            end
            if ~isempty(PDM)
                for i = 1:length(PDM)
                SEED = SEED + alpha(i+1)*fun1(PDM{i});
                end
            end
            SEED(1:N+1:end) = 0;

            B = gen_network_from_probs(SEED,desiredEdges,B);

        end

end

end