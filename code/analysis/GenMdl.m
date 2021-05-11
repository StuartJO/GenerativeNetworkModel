function [Output,E,K,P,b,C,BestFit_B,BestFit_E,BestFit_K,BestFit_C,BestFit_b,BestCorr_E,BestCorr_K,BestCorr_C,BestCorr_b] = GenMdl(A,A_dist,D,PD,Input)

% A = Adjaceny matrix to model

% A_dist = The distance matrix of A

% D = Either a single distance matrix (of some kind) for A, or a cell where
% each element contains a distance matrix for A. The model will progress
% through to the next matrix in the cell after it has added m/n edges into
% the model, where m is the total number of edges of the empirical network
% being modelled and n is the number of distance matrices

% PD = Either a single distance matrix (of some kind) for A

% Input.AddMult = 'Add' or 'Mult'. Additive or multiplicative form
% respectively. The equation forms are:
% Additive: exp(eta*D) + ( alpha1 * (T^gamma) ) +  alpha2* (PD^lambda) 
% Multiplicative: exp(eta*D) * (T^gamma) ) * (PD^lambda) 

% Input.ModelNum = A number between 1 and 14. Each corresponds to a
% different topology model:
% 1 = 'sptl'
% 2 = 'neighbors'
% 3 = 'matching'
% 4 = 'clu-avg'
% 5 = 'clu-min'
% 6 = 'clu-max'
% 7 = 'clu-diff'
% 8 = 'clu-prod'
% 9 = 'deg-avg'
% 10 = 'deg-min'
% 11 = 'deg-max'
% 12 = 'deg-diff'
% 13 = 'deg-prod'
% 14 = 'com'

% Input.Growth = 0 or 1. 0 For the static model, 1 for the growth model

% Input.DistFunc = 'powerlaw' or exponential'

% Input.GeneFunc = 'powerlaw' or exponential'

% Input.ParamRange = the range of parameters to search through as a 5*2
% matrix:
% Input.ParamRange(1,:) = eta value range
% Input.ParamRange(2,:) = gamma value range
% Input.ParamRange(3,:) = first alpha value range
% Input.ParamRange(4,:) = second alpha value range
% Input.ParamRange(5,:) = lambda alpha value range
% Note that for the multiplicative form, the alpha values will do nothing
%
% Output:
% Output.E = Energy/model fits for optimisation network
% Output.K = KS statistics of degree, clustering, betweenness and edge length for each optimisation network
% Output.P = Parameters used for each optimisation network
% Output.b = Matrix index of edges for each optimisation network
% Output.C = Spearman correlation of nodal degree (with empirical data) for each optimisation network
% Output.BestFit_B = Adjacency matrices for 100 networks produced using the best fitting parameters
% Output.BestFit_E = Energy/model fits for 100 networks produced using the best fitting parameters
% Output.BestFit_K = KS statistics for 100 networks produced using the best fitting parameters
% Output.BestFit_C = Spearman correlation of nodal degree (with empirical data) for 100 networks produced using the best fitting parameters
% Output.BestFit_b = Matrix index of edges for 100 networks produced using the best fitting parameters
% Output.BestCorr_E = Energy/model fits for 100 networks produced using the parameters producing the best correlation
% Output.BestCorr_K = KS statistics for 100 networks produced using the parameters producing the best correlation
% Output.BestCorr_C = Spearman correlation of nodal degree (with empirical data) for 100 networks produced using the parameters producing the best correlation
% Output.BestCorr_b = Matrix index of edges for 100 networks produced using the parameters producing the best correlation

% etaRange = Input.ParamRange(1,:);
% gamRange = Input.ParamRange(2,:);
% lamlim = Input.ParamRange(5,:);
% alphalim = Input.ParamRange(3,:);
% alpha2lim = Input.ParamRange(4,:);

rng('shuffle')

seed = [];

switch Input.AddMult
case 'Add'
Input.epsilon = 0;
case 'Mult'
Input.epsilon = 1e-6;
end

Input.normsum = 0;

if iscell(D) == 1 && length(D) > 1
    Steps = length(D);
    den = density_und(A); 
    m = (den/Steps):(den/Steps):den;
else
    [m,~,~] =  density_und(A); 
end

% ndraw = number of partameters to do at each optimisation level
% nlvl = number of optimisation levels
% pow = exponent of the powerlaw controlling the probability of sampling
% from the cells with the lowest energy (higher value = more likely)

Input.ndraw = 2000;
Input.nlvl = 5;
Input.pow = 2;

runparfor = 1;
%mtype = {'sptl','neighbors','matching','clu-avg','clu-min','clu-max','clu-diff','clu-prod','deg-avg','deg-min','deg-max','deg-diff','deg-prod','com'};

[E,K,P,b,C] = VoronoinLandScapeGene(Input,A,A_dist,D,PD,m,seed,runparfor);

       [~,I] = min(E); 

        bestParams = P(I,:);

        BestFit_B = cell(1,100);
        BestFit_E = zeros(1,100);
        BestFit_K = zeros(100,4);
        BestFit_C = zeros(1,100);
       parfor k = 1:100
        
        [BestFit_B{k},BestFit_b{k}] = RunGrowthModel(Input,D,PD,bestParams,m,seed);
        
           [BestFit_E(k), BestFit_K(k,:)] = Betzel_energy(A,A_dist,BestFit_B{k});
	BestFit_C(k) = corr(sum(BestFit_B{k})',sum(A)','Type','Spearman');
       end

       [~,I] = max(C); 

        bestParamsCorr = P(I,:);

        BestCorr_B = cell(1,100);
        BestCorr_E = zeros(1,100);
        BestCorr_K = zeros(100,4);
        BestCorr_C = zeros(1,100);
       parfor k = 1:100
            [BestCorr_B{k},BestCorr_b{k}] = RunGrowthModel(Input,D,PD,bestParamsCorr,m,seed);
           [BestCorr_E(k), BestCorr_K(k,:)] = Betzel_energy(A,A_dist,BestCorr_B{k});
	BestCorr_C(k) = corr(sum(BestCorr_B{k})',sum(A)','Type','Spearman');
       end
       
      
Output.E = E;
Output.K = K;
Output.P = P;
Output.b = b;
Output.C = C;
Output.BestFit_B = BestFit_B;
Output.BestFit_E = BestFit_E;
Output.BestFit_K = BestFit_K;
Output.BestFit_C = BestFit_C;
Output.BestFit_b = BestFit_b;
Output.BestCorr_E = BestCorr_E;
Output.BestCorr_K = BestCorr_K;
Output.BestCorr_C = BestCorr_C;
Output.BestCorr_b = BestCorr_b;

