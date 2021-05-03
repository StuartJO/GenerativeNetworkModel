function [EcrossValid,CcrossValid,Net] = CrossValidateModel(adjs,A_dist,D,PD,P,INSTANCES,Input)

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
addpath /projects/kg98/stuarto/BCT

for i = 1:100
%load('Rand200_GeneData.mat')
 A = double(adjs{i}>0);
 %Input.GeneFunc = 'exponential';

if Input.Growth && strcmp(Input.PD1type,'dist') == 1
    den =  density_und(A); 
    m = (den/19):(den/19):den;
else
    [m,~,~] =  density_und(A); 
end

mtype = {'sptl','neighbors','matching','clu-avg','clu-min','clu-max','clu-diff','clu-prod','deg-avg','deg-min','deg-max','deg-diff','deg-prod','com'};

for j = 1:100       

    bestEta = P(j,1);
	bestGam = P(j,2);
    bestAlpha1 = P(j,3);
    bestAlpha2 = P(j,4);    
    bestLam = P(j,5);
    
        if INSTANCES == 1

       for k = 1:INSTANCES
	   %[B{k},b_best{k}] = TopologyInput.GrowthModelAdditive(D,mtype{modelnum},bestEta,bestGam,[1 bestAlpha1 bestAlpha2],m,...
       %     1,0,seed,distfunc,'powerlaw',PD,'powerlaw',bestLam);
        
        [B,b{k}] = GrowthModel(Input.AddMult,D,mtype{Input.ModelNum},[bestEta bestLam],bestGam,[1 bestAlpha2],bestAlpha1,m,...
                1,epsilon,seed,Input.DistFunc,'powerlaw',PD,Input.GeneFunc,normsum);
        
       etemp(:,k) = Betzel_energy(A,A_dist,B);
	   ctemp(:,k) = corr(sum(B)',sum(A)','Type','Spearman');

       end            
            
        else
    

       parfor k = 1:INSTANCES
	   %[B{k},b_best{k}] = TopologyInput.GrowthModelAdditive(D,mtype{modelnum},bestEta,bestGam,[1 bestAlpha1 bestAlpha2],m,...
       %     1,0,seed,distfunc,'powerlaw',PD,'powerlaw',bestLam);
        
        [B,b{k}] = GrowthModel(Input.AddMult,D,mtype{Input.ModelNum},[bestEta bestLam],bestGam,[1 bestAlpha2],bestAlpha1,m,...
                1,epsilon,seed,Input.DistFunc,'powerlaw',PD,Input.GeneFunc,normsum);
        
       etemp(:,k) = Betzel_energy(A,A_dist,B);
	   ctemp(:,k) = corr(sum(B)',sum(A)','Type','Spearman');

       end

        end
       
       EcrossValid{j}(i,:) = etemp;
       CcrossValid{j}(i,:) = ctemp;
       Net{j,i} = b;

end
display(['Finished subject ',num2str(i)])
end