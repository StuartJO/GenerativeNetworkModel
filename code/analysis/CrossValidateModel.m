function CV_output = CrossValidateModel(adjs,A_dist,D,PD,P,INSTANCES,Input)

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

if iscell(D) == 1
    Steps = length(D);
    %den = density_und(A); 
    %m = (den/Steps):(den/Steps):den;
    [~,~,nedges] = density_und(A);
    m = round(linspace(nedges/Steps,nedges,Steps));
else
    %[m,~,~] =  density_und(A); 
    [~,~,m] =  density_und(A); 
end

mtype = {'sptl','neighbors','matching','clu-avg','clu-min','clu-max','clu-diff','clu-prod','deg-avg','deg-min','deg-max','deg-diff','deg-prod','com'};
A_nodedist = sum(A_dist.*A)./sum(A,1);
for j = 1:100       

    bestEta = P(j,1);
	bestGam = P(j,2);
    bestAlpha1 = P(j,3);
    bestAlpha2 = P(j,4);    
    bestLam = P(j,5);

        if INSTANCES == 1

       for k = 1:INSTANCES
        
        [B,b{k}] = GrowthModel(Input.AddMult,D,mtype{Input.ModelNum},[bestEta bestLam],bestGam,[1 bestAlpha2],bestAlpha1,m,...
                1,epsilon,seed,Input.DistFunc,'powerlaw',PD,Input.GeneFunc,normsum);
        
        [ksmax(k,:),KSvals(k,:),~,~,A_topo_temp,NodeVals{k}] = Betzel_energy(A,A_dist,B);
                deg(k,:) = NodeVals{k}{1};
                bet(k,:) = NodeVals{k}{2};
                clu(k,:) = NodeVals{k}{3};
                elen(k,:) = sum(A_dist.*B)./deg(k,:);
                for iii = 1:3
                ctemp(k,iii) = corr(A_topo_temp{iii},NodeVals{k}{iii},'Type','Spearman');
                end
                ctemp(k,4) = corr(A_nodedist',elen(k,:)','Type','Spearman');
                
       %etemp(:,k) = Betzel_energy(A,A_dist,B);
	   %ctemp(:,k) = corr(sum(B)',sum(A)','Type','Spearman');

       end            
          
       
        CV_output.Topography{i,j}{1} = deg;
        CV_output.Topography{i,j}{2} = bet;
        CV_output.Topography{i,j}{3} = clu;
        CV_output.Topography{i,j}{4} = elen;   
        CV_output.KS{i,j} = KSvals;
        CV_output.TopographyCorrs{i,j} = ctemp;
        CV_output.KSmax{i,j} = ksmax;

        CV_output.Net{i,j} = b;
       
        else
    

       parfor k = 1:INSTANCES
        
        [B,b{k}] = GrowthModel(Input.AddMult,D,mtype{Input.ModelNum},[bestEta bestLam],bestGam,[1 bestAlpha2],bestAlpha1,m,...
                1,epsilon,seed,Input.DistFunc,'powerlaw',PD,Input.GeneFunc,normsum);
        
        [ksmax(k,:),KSvals(k,:),~,~,A_topo_temp,NodeVals{k}] = Betzel_energy(A,A_dist,B);
                deg(k,:) = NodeVals{k}{1};
                bet(k,:) = NodeVals{k}{2};
                clu(k,:) = NodeVals{k}{3};
                elen(k,:) = sum(A_dist.*B)./deg(k,:);
                for iii = 1:4
                
                if iii == 4
                  ctemp(k,iii) = corr(A_nodedist',elen(k,:)','Type','Spearman');  
                else
                  ctemp(k,iii) = corr(A_topo_temp{k}{iii},NodeVals{k}{iii},'Type','Spearman');  
                end
                end
                
                
       end
        CV_output.Topography{i,j}{1} = deg;
        CV_output.Topography{i,j}{2} = bet;
        CV_output.Topography{i,j}{3} = clu;
        CV_output.Topography{i,j}{4} = elen;   
        CV_output.KS{i,j} = KSvals;
        CV_output.TopographyCorrs{i,j} = ctemp;
        CV_output.KSmax{i,j} = ksmax;

        CV_output.Net{i,j} = b;
        end
       


end
display(['Finished subject ',num2str(i)])
end
CV_output.Input = Input;
