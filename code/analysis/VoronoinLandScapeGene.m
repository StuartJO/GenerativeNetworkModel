function [E,K,ptsout,b,C] = VoronoinLandScapeGene(Input,A,A_dist,D,PD,density_val,InitialSEED,runparfor)

%-------------------------------------------------------------------------------
ADeg = sum(A);
numProperties = 4;

pow = Input.pow;
nlvl = Input.nlvl;
ndraw = Input.ndraw;

totalsamples = ndraw*nlvl;
ptsout = zeros(totalsamples,5);
E = zeros(totalsamples,1);
C = zeros(totalsamples,1);
K = zeros(totalsamples,numProperties);

etalim = Input.ParamRange(1,:);
gamlim = Input.ParamRange(2,:);
lamlim = Input.ParamRange(5,:);
alphalim = Input.ParamRange(3,:);
alpha2lim = Input.ParamRange(4,:);

eta = unifrnd(etalim(1),etalim(2),ndraw,1);
gam = unifrnd(gamlim(1),gamlim(2),ndraw,1);
alpha = unifrnd(alphalim(1),alphalim(2),ndraw,1);
alpha2 = unifrnd(alpha2lim(1),alpha2lim(2),ndraw,1);
lam = unifrnd(lamlim(1),lamlim(2),ndraw,1);
ptsout(1:ndraw,:) = [eta,gam,alpha,alpha2,lam];
bounds = [etalim;gamlim;alphalim;alpha2lim;lamlim];
%-------------------------------------------------------------------------------
powvals = linspace(0,pow,nlvl);

TopoTypes = {'sptl','neighbors','matching','clu-avg','clu-min','clu-max','clu-diff','clu-prod','deg-avg','deg-min','deg-max','deg-diff','deg-prod','com'};

TopoType = TopoTypes{Input.ModelNum};

%b = [];
b = cell(1,totalsamples);
for ilvl = 1:nlvl
    fprintf('level %i of %i\n',ilvl,nlvl);

    pow = powvals(ilvl);

    % If the first level, use the randomly selected points. If > the first
    % level, run the voronoi tessellation to sample new parameters from the
    % cells with the lowest energy/bets model fits
    
    if ilvl==1
        ptsnew = ptsout(1:ndraw,:);
    else
        ind = 1:ndraw*(ilvl - 1);
        if strcmp(TopoType,'sptl')
		if isempty(PD)
            ptsnew = [fcn_voronoi_select_sptl(ptsout(ind,1),E(ind),ndraw,etalim,pow),gam,alpha,alpha2,lam];
		else
            ptsnew = fcn_voronoi_selectn(ptsout(ind,:),E(ind),ndraw,bounds,pow);
		end
        else
            
            ptsnew = fcn_voronoi_selectn(ptsout(ind,:),E(ind),ndraw,bounds,pow);
        end
    end

    indnew = (1 + (ilvl - 1)*ndraw):(ilvl*ndraw);
    
    % Run in a parfor loop or not
    
    if runparfor
       
        Epar = zeros(1,ndraw);
        Kpar = zeros(ndraw,4);
        Cpar = zeros(1,ndraw);
        btemp = cell(1,ndraw);
                
        parfor i = 1:ndraw     
                [B,btemp{i}] = RunGrowthModel(Input,D,PD,ptsnew(i,:),density_val,InitialSEED);
            [Epar(i),Kpar(i,:)] = Betzel_energy(A,A_dist,B);
            Cpar(i) = corr(sum(B)',ADeg','Type','Spearman');
        end
        
        E(indnew) = Epar;
        K(indnew,:) = Kpar;
        C(indnew) = Cpar;
        b(indnew) = btemp;
        
    else

    for i = 1:ndraw
        indHere = indnew(i);
                [B,b{i}]  = RunGrowthModel(Input,D,PD,ptsnew(i,:),density_val,InitialSEED);
        [E(indHere),K(indHere,:)] = Betzel_energy(A,A_dist,B);
        C(indHere) = corr(sum(B)',ADeg','Type','Spearman');

    end
    
    end

    ptsout(indnew,:) = ptsnew;
end
