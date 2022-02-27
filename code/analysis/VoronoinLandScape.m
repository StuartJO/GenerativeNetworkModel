function [maxKS,KS,P,b,C] = VoronoinLandScape(A,A_dist,PD1,PD2,m,Input)

% This function will perform optimisation using a Voronoi tessellation
% approach, where from an initial set of random points in parameter space,
% Voronoi tessellation will be performed to divide up the parameter space
% into cells. Each cell has an associated maxKS value (based on the
% parameter value used to draw that cell). The algorithm proceeds by
% preferentially sampling parameter values from cells with the smallest 
% maxKS, and then performing Voronoi tessellation again. This repeats a
% number of times.

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
C = zeros(totalsamples,1);
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
                [maxKSpar(i),KSpar(i,:)] = Betzel_energy(A,A_dist,B);
                Cpar(i) = corr(sum(B)',ADeg','Type','Spearman');
        end
        
        maxKS(indnew) = maxKSpar;
        KS(indnew,:) = KSpar;
        C(indnew) = Cpar;
        b(indnew) = btemp;
        
    else

        for i = 1:ndraw
            indHere = indnew(i);
            [B,b{indHere}] = GrowthModel(PD1,PD2,ptsnew(i,:),m,Input);
            [maxKS(indHere),KS(indHere,:)] = Betzel_energy(A,A_dist,B);
            C(indHere) = corr(sum(B)',ADeg','Type','Spearman');
        end
    
    end

    P(indnew,:) = ptsnew;
end
