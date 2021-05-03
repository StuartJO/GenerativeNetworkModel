function [Output,E,K,P,b,C,B,Ebest,Kbest,Cbest,b_best,EbestCorr,KbestCorr,CbestCorr,b_bestCorr] = GenMdl(A,A_dist,D,PD,Input)

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
% Input.ParamRange(4,:) = lambda alpha value range
% Note that for the multiplicative form, the alpha values will do nothing

etaRange = Input.ParamRange(1,:);
gamRange = Input.ParamRange(2,:);
lamlim = Input.ParamRange(5,:);
alphalim = Input.ParamRange(3,:);
alpha2lim = Input.ParamRange(4,:);

rng('shuffle')

seed = [];
 %end
switch Input.AddMult
case 'Add'
epsilon = 0;
case 'Mult'
epsilon = 1e-6;
end

normsum = 0;


if iscell(D) == 1 && length(D) > 1
    Steps = length(D);
    den = density_und(A); 
    m = (den/Steps):(den/Steps):den;
else
    [m,~,~] =  density_und(A); 
end

ndraw = 2000;

mtype = {'sptl','neighbors','matching','clu-avg','clu-min','clu-max','clu-diff','clu-prod','deg-avg','deg-min','deg-max','deg-diff','deg-prod','com'};

[E,K,P,b,C] = VoronoinLandScapeGene(Input.AddMult,A,A_dist,D,mtype{Input.ModelNum},m,seed,Input.DistFunc,'powerlaw',etaRange,gamRange,alphalim,2,5,ndraw,1,PD,lamlim,alpha2lim,Input.GeneFunc,normsum);


       [~,I] = min(E); 
       

    bestEta = P(I,1);
	bestGam = P(I,2);
    bestAlpha1 = P(I,3);
    bestAlpha2 = P(I,4);    
    bestLam = P(I,5);

        B = cell(1,100);
        Ebest = zeros(1,100);
        Kbest = zeros(100,4);
        Cbest = zeros(1,100);
       parfor k = 1:100
        
        [B{k},b_best{k}] = GrowthModel(Input.AddMult,D,mtype{Input.ModelNum},[bestEta bestLam],bestGam,[1 bestAlpha2],bestAlpha1,m,...
                1,epsilon,seed,Input.DistFunc,'powerlaw',PD,Input.GeneFunc,normsum);
        
           [Ebest(k), Kbest(k,:)] = Betzel_energy(A,A_dist,B{k});
	Cbest(k) = corr(sum(B{k})',sum(A)','Type','Spearman');
       end



       [~,I] = max(C); 

    bestEtaCorr = P(I,1);
	bestGamCorr = P(I,2);
    bestAlpha1Corr = P(I,3);
    bestAlpha2Corr = P(I,4);
    bestLamCorr = P(I,5);

       modeltype = mtype{Input.ModelNum};
        BCorr = cell(1,100);
        EbestCorr = zeros(1,100);
        KbestCorr = zeros(100,4);
        CbestCorr = zeros(1,100);
       parfor k = 1:100
         %  [BCorr{k},b_bestCorr{k}] = TopologyInput.GrowthModelAdditive(D,mtype{modelnum},bestEtaCorr,bestGamCorr,[1 bestAlphaCorr],m,...
           % 1,0,seed,distfunc,'powerlaw');
                [BCorr{k},b_bestCorr{k}]  = GrowthModel(Input.AddMult,D,mtype{Input.ModelNum},[bestEtaCorr bestLamCorr],bestGamCorr,[1 bestAlpha2Corr],bestAlpha1Corr,m,...
                1,epsilon,seed,Input.DistFunc,'powerlaw',PD,Input.GeneFunc,normsum);
           [EbestCorr(k), KbestCorr(k,:)] = Betzel_energy(A,A_dist,BCorr{k});
	CbestCorr(k) = corr(sum(BCorr{k})',sum(A)','Type','Spearman');
       end
       
      
Output.E = E;
Output.K = K;
Output.P = P;
Output.b = b;
Output.C = C;
Output.B = B;
Output.Ebest = Ebest;
Output.Kbest = Kbest;
Output.Cbest = Cbest;
Output.b_best = b_best;
Output.EbestCorr = EbestCorr;
Output.KbestCorr = KbestCorr;
Output.CbestCorr = CbestCorr;
Output.b_bestCorr = b_bestCorr;

