function MdlOutput = GenMdl(A,A_dist,PD1,PD2,Input)

% This function runs the generative model optimisation. 
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
%
% Outputs:
% 
% MdlOutput = a structure with the following fields:
%   MdlOutput.maxKS = N*1 array of maxKS values for each network generated
%       during optimisation
%   MdlOutput.DegCorr = N*1 array of values for the correlation with  
%       empirical degree each network generated during optimisation  
%   MdlOutput.KS = N*4 matrix of KS values for each network generated
%       during optimisation. Each column corresponds to the following 
%       measure: 1 = degree; 2 = clustering; 3 = betweenness; 4 = mean edge 
%       length
%   MdlOutput.P = N*5 matrix of parameter values for each network generated
%       during optimisation
%   MdlOutput.b = 1*N cell array, where each cell contains an edge index
%       list for each network generated during optimisation
% 
%   MdlOutput.optim_maxKS = for networks generated using the parameters
%       of the lowest maxKS value, a 1*100 array of maxKS values
%   MdlOutput.optim_KS = for networks generated using the parameters
%       of the lowest maxKS value, an 100*4 matrix of KS values
%   MdlOutput.optim_b = a 1*100 cell array where each cell contains an edge
%       index list for each network generated during using the parameters 
%       of the lowest maxKS value       
%   MdlOutput.optim_DegCorr = for networks generated using the parameters
%       of the lowest maxKS value, a 1*100 array of values for the 
%       correlation with empirical degree      
%
%   MdlOutput.bestDegCorr_maxKS = for networks generated using the parameters
%       of the largest degree correlation, a 1*100 array of maxKS values
%   MdlOutput.bestDegCorr_KS = for networks generated using the parameters
%       of the largest degree correlation, an 100*4 matrix of KS values
%   MdlOutput.bestDegCorr_b = a 1*100 cell array where each cell contains an edge
%       index list for each network generated during using the parameters 
%       of the largest degree correlation    
%   MdlOutput.bestDegCorr_DegCorr = for networks generated using the parameters
%       of the largest degree correlation, a 1*100 array of values for the 
%       correlation with empirical degree     
%
% You may wonder why the model has a different form to what we report in
% the paper e.g., (D^eta)+a(T^gam) or (D^eta)+a(PC^gam) etc. Simply,
% because a) we never run a model with distance + topology + gene, b) we
% wanted to avoid discussing too many parameters because that would get
% confusing quick and readers would find it difficult to keep track of, and
% c) topology and PC require seperate configuration of their parameters the
% way the model is coded. So while the parameters in the topological and PC
% models perform the same function theoretically, practically they are
% different. This however does mean you can run a model with all 5 
% parameters...(good luck getting the optimisation to work!)
%
% Another question you may have is why the ParamRange variable is ordered in
% this way i.e., two exponents, then the two alpha values, then an
% exponent. The simple answer is this just follows the order in which the
% function developed over time. eta and gamma are the OG parameters so they
% come first. We tried a simple model using only alpha for the topology
% term so that came third. Then we added the PD2 term, which we first tried 
% with just an alpha value (4th), then with its own exponent as well (5th).

rng('shuffle')

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

% Basically this assumes you want the same number of edges formed for each
% of the input distance matrices.

if iscell(PD1) == 1 && length(PD1) > 1
    Steps = length(PD1);
    [~,~,nedges] = density_und(A);
    m = round(linspace(nedges/Steps,nedges,Steps));
else
    [~,~,m] =  density_und(A); 
end

[maxKS,KS,P,b,DegCorr] = VoronoinLandScape(A,A_dist,PD1,PD2,m,Input);

[~,I] = min(maxKS);     
P_optim = P(I,:);
optim_b = cell(1,100);

optim_maxKS = zeros(1,100);
optim_KS = zeros(100,4);
optim_DegCorr = zeros(1,100);

if Input.useParfor

    parfor k = 1:100
        [B,optim_b{k}] = GrowthModel(PD1,PD2,P_optim,m,Input);
        [optim_maxKS(k), optim_KS(k,:)] = calc_maxKS(A,A_dist,B);
        optim_DegCorr(k) = corr(sum(B)',sum(A)','Type','Spearman');
    end   
else
    for k = 1:100
        [B,optim_b{k}] = GrowthModel(PD1,PD2,P_optim,m,Input);
        [optim_maxKS(k), optim_KS(k,:)] = calc_maxKS(A,A_dist,B);
        optim_DegCorr(k) = corr(sum(B)',sum(A)','Type','Spearman');
    end     
end

[~,I] = max(DegCorr);     
bestDegCorr_P = P(I,:);
bestDegCorr_b = cell(1,100);

bestDegCorr_maxKS = zeros(1,100);
bestDegCorr_KS = zeros(100,4);
bestDegCorr_DegCorr = zeros(1,100);

if Input.useParfor

    parfor k = 1:100
        [B,bestDegCorr_b{k}] = GrowthModel(PD1,PD2,bestDegCorr_P,m,Input);
        [bestDegCorr_maxKS(k), bestDegCorr_KS(k,:)] = calc_maxKS(A,A_dist,B);
        bestDegCorr_DegCorr(k) = corr(sum(B)',sum(A)','Type','Spearman');
    end 

else
    
    for k = 1:100
        [B,bestDegCorr_b{k}] = GrowthModel(PD1,PD2,bestDegCorr_P,m,Input);
        [bestDegCorr_maxKS(k), bestDegCorr_KS(k,:)] = calc_maxKS(A,A_dist,B);
        bestDegCorr_DegCorr(k) = corr(sum(B)',sum(A)','Type','Spearman');
    end  

end
       
MdlOutput.maxKS = maxKS;
MdlOutput.DegCorr = DegCorr;
MdlOutput.KS = KS;
MdlOutput.P = P;
MdlOutput.b = b;

MdlOutput.optim_maxKS = optim_maxKS;
MdlOutput.optim_KS = optim_KS;
MdlOutput.optim_b = optim_b;
MdlOutput.optim_DegCorr = optim_DegCorr;

MdlOutput.bestDegCorr_maxKS = bestDegCorr_maxKS;
MdlOutput.bestDegCorr_KS = bestDegCorr_KS;
MdlOutput.bestDegCorr_b = bestDegCorr_b;
MdlOutput.bestDegCorr_DegCorr = bestDegCorr_DegCorr;

% Save the input configurations to output. Helps to keep track of what was
% done
MdlOutput.Input = Input;
MdlOutput.Input.Nnodes = length(A);

