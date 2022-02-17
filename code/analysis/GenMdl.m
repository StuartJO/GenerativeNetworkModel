function Output = GenMdl(A,A_dist,PD1,PD2,Input)

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
% Teachnically, this function is a wrapper for another function
% (GrowthModel) which itself is a wrapper for functions which actually run
% the generative model. Turtles all the way down.
%
% A number of inputs are required:
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
%   topology model:
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
%  ParamRange = A 5x2 matrix where the first column gives the lower, and
%  the second the upper, bounds for the following parameters (see the basic
%  form outlined at the top):
%       ParamRange(1,:) = eta
%       ParamRange(2,:) = gam
%       ParamRange(3,:) = a1
%       ParamRange(4,:) = a2
%       ParamRange(5,:) = lam
%
% The following are optional specifications for Input (will be set to
% defaults):
%
% PD1_Func = 'power-law' or 'exponential'. Controls the interaction between
% PD1 and eta
%
% T_Func = 'power-law' or 'exponential'. Controls the interaction between T
% and eta
%
% PD2_Func = 'power-law' or 'exponential'. Controls the interaction between
% PD2 and lam
%
% See the bottom 
%
% You may wonder why there are 5 possible parameters but we only ever talk
% about a max of three in the paper. Very simply, topology and PC require
% seperate configuration of their parameters, and it was far easier to keep
% track of what is what if all parameter variables were set up the same. SO
% while functionally the parameters in the topological and PC models are
% the same, practically they appear different. This however does mean you
% can do a model with all 5 parameters...

% Another question you may have is why the ParamRange variable is ordered in
% this way i.e., two exponents, then the two alpha values, then an
% exponent. The simple answer is this just follows the order in which the
% function developed over time. eta and gamma are the OG parameters so they
% come first. We tried a simple model using only alpha for the topology
% term so that came third. Then we added the PD term, which we first tried 
% with just an alpha value (4th), then with its own exponent as well (5th).
% In the paper 
etaRange = Input.ParamRange(1,:);
gamRange = Input.ParamRange(2,:);
lamRange = Input.ParamRange(5,:);
a1Range = Input.ParamRange(3,:);
a2Range = Input.ParamRange(4,:);

rng('shuffle')

 %end
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

TopoTypes = {'sptl','neighbors','matching','clu-avg','clu-min','clu-max','clu-diff','clu-prod','deg-avg','deg-min','deg-max','deg-diff','deg-prod','com'};
Input.TopoType = TopoTypes{Input.ModelNum};

% Basically this assumes you want the same number of edges formed for each
% of the input distance matrices.

if iscell(PD1) == 1 && length(PD1) > 1
    Steps = length(PD1);
    [~,~,nedges] = density_und(A);
    m = round(linspace(nedges/Steps,nedges,Steps));
else
    [~,~,m] =  density_und(A); 
end

%Input.ndraw = 10;
% Input.ndraw = 2000;
% Input.pow = 2;
% Input.nlvl = 5;
% Input.useParFor = 1;

[maxKS,KS,P,b] = VoronoinLandScape(A,A_dist,PD1,PD2,m,Input);

    [~,I] = min(maxKS);     
    P_optim = P(I,:);
    b_optim = cell(1,100);

        maxKS_optim = zeros(1,100);
        KS_optim = zeros(100,4);
        %C_optim = zeros(1,100);
       parfor k = 1:100
        [B,b_optim{k}] = GrowthModel(PD1,PD2,P_optim,m,Input);
           [maxKS_optim(k), KS_optim(k,:)] = Betzel_energy(A,A_dist,B);
            %C_optim(k) = corr(sum(B)',sum(A)','Type','Spearman');
       end      
      
Output.maxKS = maxKS;
Output.KS = KS;
Output.P = P;
Output.b = b;
Output.maxKS_optim = maxKS_optim;
Output.KS_optim = KS_optim;
Output.b_optim = b_optim;

% Save the input configurations to output. Helps to keep track of what was
% done
Output.Input = Input;

