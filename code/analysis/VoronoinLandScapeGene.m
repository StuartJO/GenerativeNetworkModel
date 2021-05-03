function [E,K,ptsout,b,C] = VoronoinLandScapeGene(AddMult,A,A_dist,D,modeltype,den_levels,InitialSEED,DistFunc,TopoFunc,etalim,gamlim,alphalim,pow,nlvl,ndraw,runparfor,O,lamlim,alpha2lim,GeneFunc,normsum)

% DistFunc = 'powerlaw';
% TopoFunc = 'powerlaw';
if nargin < 20
    GeneFunc = 'powerlaw';
end

if nargin < 21
    normsum = 0;
end

switch AddMult
case 'Add'
epsilon = 0;
case 'Mult'

epsilon = 1e-5;
end
%-------------------------------------------------------------------------------
ADeg = sum(A);
numProperties = 4;

totalsamples = ndraw*nlvl;
ptsout = zeros(totalsamples,5);
E = zeros(totalsamples,1);
C = zeros(totalsamples,1);
K = zeros(totalsamples,numProperties);

eta = unifrnd(etalim(1),etalim(2),ndraw,1);
gam = unifrnd(gamlim(1),gamlim(2),ndraw,1);
alpha = unifrnd(alphalim(1),alphalim(2),ndraw,1);
alpha2 = unifrnd(alpha2lim(1),alpha2lim(2),ndraw,1);
lam = unifrnd(lamlim(1),lamlim(2),ndraw,1);
ptsout(1:ndraw,:) = [eta,gam,alpha,alpha2,lam];
bounds = [etalim;gamlim;alphalim;alpha2lim;lamlim];
%-------------------------------------------------------------------------------
powvals = linspace(0,pow,nlvl);
bestE = 1;
%b = [];
b = cell(1,totalsamples);
for ilvl = 1:nlvl
    fprintf('level %i of %i\n',ilvl,nlvl);

    pow = powvals(ilvl);

    if ilvl==1
        ptsnew = ptsout(1:ndraw,:);
    else
        ind = 1:ndraw*(ilvl - 1);
        if strcmp(modeltype,'sptl')
            %ptsnew = [fcn_voronoi_select_sptl(ptsout(ind,1),E(ind),ndraw,etalim,pow),gam,alpha,alpha2,lam];
		if isempty(O)
		ptsnew = [fcn_voronoi_select_sptl(ptsout(ind,1),E(ind),ndraw,etalim,pow),gam,alpha,alpha2,lam];
		else
            ptsnew = fcn_voronoi_selectn(ptsout(ind,:),E(ind),ndraw,bounds,pow);
		end
        else
            
            ptsnew = fcn_voronoi_selectn(ptsout(ind,:),E(ind),ndraw,bounds,pow);
        end
    end

    indnew = (1 + (ilvl - 1)*ndraw):(ilvl*ndraw);
    
    if runparfor
       
        Epar = zeros(1,ndraw);
        Kpar = zeros(ndraw,4);
        Cpar = zeros(1,ndraw);
        btemp = cell(1,ndraw);
        
        etaptsnew = ptsnew(:,1);
        gamptsnew = ptsnew(:,2);
        alphaptsnew = ptsnew(:,3);
        alpha2ptsnew = ptsnew(:,4);
        lamptsnew = ptsnew(:,5);
        
        parfor i = 1:ndraw          
%             [B,btemp{i}] = GrowthModel(D,modeltype,[etaptsnew(i) lamptsnew(i)],gamptsnew(i),[1 alpha2ptsnew(i)],alphaptsnew(i),den_levels,...
%                 1,0,InitialSEED,DistFunc,TopoFunc,O);
                [B,btemp{i}]  = GrowthModel(AddMult,D,modeltype,[etaptsnew(i) lamptsnew(i)],gamptsnew(i),[1 alpha2ptsnew(i)],alphaptsnew(i),den_levels,...
                1,epsilon,InitialSEED,DistFunc,TopoFunc,O,GeneFunc,normsum);
            [Epar(i),Kpar(i,:)] = Betzel_energy(A,A_dist,B);
            Cpar(i) = corr(sum(B)',ADeg','Type','Spearman');
        end
        
        E(indnew) = Epar;
        K(indnew,:) = Kpar;
        C(indnew) = Cpar;
        b(indnew) = btemp;
%        [~,I] = min(Epar);
%         if min(Epar) < bestE
%             b = btemp{I};
%             bestE = min(Epar);
%         end
        
    else
        etaptsnew = ptsnew(:,1);
        gamptsnew = ptsnew(:,2);
        alphaptsnew = ptsnew(:,3);
        alpha2ptsnew = ptsnew(:,4);
        lamptsnew = ptsnew(:,5);
    for i = 1:ndraw
        indHere = indnew(i);
%         [B,btemp] = GrowthModel(D,modeltype,[etaptsnew(i) lamptsnew(i)],gamptsnew(i),[1 alpha2ptsnew(i)],alphaptsnew(i),den_levels,...
%                 1,0,InitialSEED,DistFunc,TopoFunc,O);
                [B,btemp]  = GrowthModel(AddMult,D,modeltype,[etaptsnew(i) lamptsnew(i)],gamptsnew(i),[1 alpha2ptsnew(i)],alphaptsnew(i),den_levels,...
                1,epsilon,InitialSEED,DistFunc,TopoFunc,O,GeneFunc,normsum);
        [E(indHere),K(indHere,:)] = Betzel_energy(A,A_dist,B);
        C(indHere) = corr(sum(B)',ADeg','Type','Spearman');
        if E(indHere) < bestE
            b = btemp;
            bestE = E(indHere);
        end

    end
    
    end

    ptsout(indnew,:) = ptsnew;
end
