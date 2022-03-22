function CV_output = CrossValidateModel(adjs,A_dist,PD1,PD2,P,INSTANCES,Input)

% This function will for a set of networks, run each against a set of parameters
% to generate a network, and calculate KS, max(KS), and topography 
% information off of this
% Inputs:
% adjs = the adjacency matrices to model
%
% A_dist = The distances between nodes (e.g., Euclidean distances, fibre
% distance) for adjs. This distance is used to compute how similar the model
% is to the input adjacency matrix. A single A_dist can be provided or a
% cell array can be provided, where each cell correspond to the distances
% for the corresponding network in adjs
%
% PD1 = Either a matrix indicating pairwise similarity/distances between nodes, or a cell where each
% element contains such a distance matrix, to be used for the modelling itself.
% The "growth" model is specified by the latter of these options. If there
% is an individual PD1 for each network in adjs, PD1 can be setup as a cell
% array (much like A_dist). If A_dist is a cell, then PD1 must also be
% (even if the same PD1 matrix is applied to all subjects)
%
% PD2 = A second matrix indicating pairwise similarity/distances between nodes in
% A (unlike with PD1 a cell cannot be used as an input here). This could be 
% correlated gene expression, similarity in histology etc. As with
% A_dist/PD1, can be supplied as a cell for individualised matrices
%
% INSTANCES = how many repeats of each network and parameter combination to
% run
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
% Output:
% 
% CV_output = a structure with the following fields
%   CV_output.KS = the KS values for all networks produced during crossvalidation
%
%   CV_output.TopographyCorrs = the correlation between the model network and 
%       empirical topological data for all networks produced during
%       crossvalidation
%
%   CV_output.Fcv = the FCV value for all networks produced during crossvalidation
%       
%   CV_output.P = the parameters used during crossvalidation
%
%   CV_output.Topography{t} = the topology values for measure 't' (1 = degree
%       ; 2 = clustering; 3 = betweenness; 4 = mean edge length) for all 
%       networks produced during crossvalidation for model number 'mdl'.
%       Rows correspond to each network, columns correspond to nodes.
%
%   CV_output.Input = the initial set of input settings


% Need to shuffle if running on a cluster
rng('shuffle')

TopoTypes = {'sptl','neighbors','matching','clu-avg','clu-min','clu-max','clu-diff','clu-prod','deg-avg','deg-min','deg-max','deg-diff','deg-prod','com'};
Input.TopoType = TopoTypes{Input.ModelNum};

useParfor = Input.useParfor;

Nsubs = length(adjs);
Nparams = size(P,1);

if iscell(A_dist)
   A_dists = A_dist; 
   PD1s = PD1;
end

if iscell(PD2)
   PD2s = PD2; 
end

for i = 1:Nsubs
    
A = double(adjs{i}>0);

if iscell(PD1) == 1
    Steps = length(PD1);
    [~,~,nedges] = density_und(A);
    m = round(linspace(nedges/Steps,nedges,Steps));
else
    [~,~,m] = density_und(A); 
end

if iscell(A_dist)
   A_dist = A_dists{i}; 
   PD1 = PD1s{i};
end

if iscell(PD2)
   PD2 = PD2s{i}; 
end

A_nodedist = sum(A_dist.*A)./sum(A,1);
NNodes = length(A);

for j = 1:Nparams       
    p = P(j,:);
    
        maxks = zeros(INSTANCES,1);
        KSvals = zeros(INSTANCES,4);
        deg = zeros(INSTANCES,NNodes);
        clu = zeros(INSTANCES,NNodes);
        bet = zeros(INSTANCES,NNodes);
        elen = zeros(INSTANCES,NNodes);
        ctemp = zeros(INSTANCES,4);
        b = cell(1,INSTANCES);
        
        % Don't bother with parfor if only making one iteration
        if ~useParfor || INSTANCES == 1

         for k = 1:INSTANCES            
        
            [B,b{k}] = GrowthModel(PD1,PD2,p,m,Input);

            [maxks(k,:),KSvals(k,:),~,~,A_topo_temp,NodeVals] = calc_maxKS(A,A_dist,B);
            deg(k,:) = NodeVals{1};
            clu(k,:) = NodeVals{2};
            bet(k,:) = NodeVals{3};
            elen(k,:) = sum(A_dist.*B)./deg(k,:);
            for iii = 1:3
                ctemp(1,iii) = corr(A_topo_temp{iii},NodeVals{iii},'Type','Spearman');
            end
            ctemp(k,4) = corr(A_nodedist',elen(k,:)','Type','Spearman');
                       
        CV_output.Topography{i,j}{1} = deg;
        CV_output.Topography{i,j}{2} = clu;
        CV_output.Topography{i,j}{3} = bet;
        CV_output.Topography{i,j}{4} = elen;   
        CV_output.KS{i,j} = KSvals;
        CV_output.TopographyCorrs{i,j} = ctemp;
        CV_output.KSmax{i,j} = maxks;

        CV_output.net{i,j} = b;
        
         end
       
        else
    
       parfor k = 1:INSTANCES
        
        [B,b{k}] = GrowthModel(PD1,PD2,p,m,Input);
        
        [maxks(k,:),KSvals(k,:),~,~,A_topo_temp,NodeVals] = calc_maxKS(A,A_dist,B);
                deg(k,:) = NodeVals{1};
                clu(k,:) = NodeVals{2};
                bet(k,:) = NodeVals{3};
                elen(k,:) = sum(A_dist.*B)./deg(k,:);
                for iii = 1:4     
                    if iii == 4
                      ctemp(k,iii) = corr(A_nodedist',elen(k,:)','Type','Spearman');  
                    else
                      ctemp(k,iii) = corr(A_topo_temp{iii},NodeVals{iii},'Type','Spearman');  
                    end
                end
       
       end
        CV_output.Topography{i,j}{1} = deg;
        CV_output.Topography{i,j}{2} = clu;
        CV_output.Topography{i,j}{3} = bet;
        CV_output.Topography{i,j}{4} = elen;   
        CV_output.KS{i,j} = KSvals;
        CV_output.TopographyCorrs{i,j} = ctemp;
        CV_output.maxKS{i,j} = maxks;

        CV_output.net{i,j} = b;
        end
       
end
display(['Finished subject ',num2str(i)])
end
CV_output.Input = Input;
CV_output.Input.NNodes = NNodes;
