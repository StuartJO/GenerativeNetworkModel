function CV_output = CrossValidateModel(adjs,A_dist,PD1,PD2,P,INSTANCES,Input)

% This function will for a set of networks, run each against a set of parameters
% to generate a network, and calculate KS, max(KS), and topography 
% information off of this

rng('shuffle')

TopoTypes = {'sptl','neighbors','matching','clu-avg','clu-min','clu-max','clu-diff','clu-prod','deg-avg','deg-min','deg-max','deg-diff','deg-prod','com'};
Input.TopoType = TopoTypes{Input.ModelNum};

useParfor = Input.useParfor;

Nsubs = length(adjs);
Nparams = size(P,1);

for i = 1:Nsubs
    
A = double(adjs{i}>0);

if iscell(PD1) == 1
    Steps = length(PD1);
    [~,~,nedges] = density_und(A);
    m = round(linspace(nedges/Steps,nedges,Steps));
else
    [~,~,m] = density_und(A); 
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
        
        if ~useParfor || INSTANCES == 1

         for k = 1:INSTANCES            
        
            [B,b{k}] = GrowthModel(PD1,PD2,p,m,Input);

            [maxks(k,:),KSvals(k,:),~,~,A_topo_temp,NodeVals] = Betzel_energy(A,A_dist,B);
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
        
        [maxks(k,:),KSvals(k,:),~,~,A_topo_temp,NodeVals] = Betzel_energy(A,A_dist,B);
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
