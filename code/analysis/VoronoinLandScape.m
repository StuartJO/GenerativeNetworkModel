function [maxKS,KS,P,b,DegCorr] = VoronoinLandScape(A,A_dist,PD1,PD2,m,Input)

% This function will perform optimisation using a Voronoi tessellation
% approach, where from an initial set of random points in parameter space,
% Voronoi tessellation will be performed to divide up the parameter space
% into cells. Each cell has an associated maxKS value (based on the
% parameter value used to draw that cell). The algorithm proceeds by
% preferentially sampling parameter values from cells with the smallest 
% maxKS, and then performing Voronoi tessellation again. This repeats a
% number of times.
%
% Inputs:
%
% A = The adjacency matrix to model
%
% A_dist = The distances between nodes (e.g., Euclidean distances, fibre
% distance) for A. This distance is used to compute how similar the model
% is to the input adjacency matrix.
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
%  pow = the severity/exponent for the opotimisation. Higher values will
%  make it more likely to sample from cells that have the smallest maxKS
%
%  nlvl = number of steps in the optimisation
%
%  ndraw = number of repetitions/samples per step
%
%   The following are optional specifications for Input (will be set to
%   defaults if not set):
%
%   PD1Func = 'power-law' or 'exponential'. Controls the interaction between
%   PD1 and eta. Defaults to 'exponential'
%
%   TopoFunc = 'power-law' or 'exponential'. Controls the interaction between T
%   and eta. Defaults to 'power-law'
%
%   PD2Func = 'power-law' or 'exponential'. Controls the interaction between
%   PD2 and lam. Defaults to 'power-law'
%
%   epsilon = an amount to add to each edges topology value in the model (to
%   ensure each edge doesn't become undefinied). Defaults to 0 when 
%   Input.AddMult = 'Add' or 1e-6 when Input.AddMult = 'Mult'
%
%   seed = a seed network. If none is desired set to [] (default)
%
%   useParfor = set to 1 to use parfor loops where possible (0 by default)
%
%   normsum = set to 1 to normalise each term by its sum (set to 0 by
%   default, which normalises by the max)


if ~isfield(Input,'epsilon')
    switch Input.AddMult
    case 'Add'
    Input.epsilon = 0;
    case 'Mult'
    Input.epsilon = 1e-6;
    end
end

if ~isfield(Input,'normsum')
    Input.normsum = 0;  
end

if ~isfield(Input,'seed')
    Input.seed = [];  
end

if ~isfield(Input,'useParfor')
 Input.useParfor = 0;  
end

TopoTypes = {'sptl','neighbors','matching','clu-avg','clu-min','clu-max','clu-diff','clu-prod','deg-avg','deg-min','deg-max','deg-diff','deg-prod','com'};
if ~isfield(Input,'TopoTypes') 
    Input.TopoType = TopoTypes{Input.ModelNum};
else
   if ~strcmp(TopoTypes{Input.ModelNum},Input.TopoType) 
      error('Input.ModelNum and Input.TopoType do not match!')
   end
end

ADeg = sum(A);
numProperties = 4;

pow = Input.pow;
nlvl = Input.nlvl;
ndraw = Input.ndraw;

useParfor = Input.useParfor;

TopoType = Input.TopoType;

totalsamples = ndraw*nlvl;

% Initalise output variables
P = zeros(totalsamples,5);
maxKS = zeros(totalsamples,1);
DegCorr = zeros(totalsamples,1);
KS = zeros(totalsamples,numProperties);
powvals = linspace(0,pow,nlvl);
b = cell(1,totalsamples);

% Pull out the different parameters and sample at random from their ranges
% to get a set of points to start the optimisation from
etaRange = Input.ParamRange(1,:);
gamRange = Input.ParamRange(2,:);
lamRange = Input.ParamRange(5,:);
a1Range = Input.ParamRange(3,:);
a2Range = Input.ParamRange(4,:);

eta = unifrnd(etaRange(1),etaRange(2),ndraw,1);
gam = unifrnd(gamRange(1),gamRange(2),ndraw,1);
a1 = unifrnd(a1Range(1),a1Range(2),ndraw,1);
a2 = unifrnd(a2Range(1),a2Range(2),ndraw,1);
lam = unifrnd(lamRange(1),lamRange(2),ndraw,1);
P(1:ndraw,:) = [eta,gam,a1,a2,lam];
bounds = [etaRange;gamRange;a1Range;a2Range;lamRange];

for ilvl = 1:nlvl
    fprintf('level %i of %i\n',ilvl,nlvl);

    pow = powvals(ilvl);

    if ilvl==1
        ptsnew = P(1:ndraw,:);
    else
        ind = 1:ndraw*(ilvl - 1);
        if strcmp(TopoType,'sptl')
		if isempty(PD2)
            ptsnew = [fcn_voronoi_select_sptl(P(ind,1),maxKS(ind),ndraw,etaRange,pow),gam,a1,a2,lam];
		else
            ptsnew = fcn_voronoi_selectn(P(ind,:),maxKS(ind),ndraw,bounds,pow);
		end
        else
            ptsnew = fcn_voronoi_selectn(P(ind,:),maxKS(ind),ndraw,bounds,pow);
        end
    end

    indnew = (1 + (ilvl - 1)*ndraw):(ilvl*ndraw);
    
    if useParfor
       
        maxKSpar = zeros(1,ndraw);
        KSpar = zeros(ndraw,4);
        Cpar = zeros(1,ndraw);
        btemp = cell(1,ndraw);
                
        parfor i = 1:ndraw          
                [B,btemp{i}] = GrowthModel(PD1,PD2,ptsnew(i,:),m,Input);
                [maxKSpar(i),KSpar(i,:)] = calc_maxKS(A,A_dist,B);
                Cpar(i) = corr(sum(B)',ADeg','Type','Spearman');
        end
        
        maxKS(indnew) = maxKSpar;
        KS(indnew,:) = KSpar;
        DegCorr(indnew) = Cpar;
        b(indnew) = btemp;
        
    else

        for i = 1:ndraw
            indHere = indnew(i);
            [B,b{indHere}] = GrowthModel(PD1,PD2,ptsnew(i,:),m,Input);
            [maxKS(indHere),KS(indHere,:)] = calc_maxKS(A,A_dist,B);
            DegCorr(indHere) = corr(sum(B)',ADeg','Type','Spearman');
        end
    
    end

    P(indnew,:) = ptsnew;
end
