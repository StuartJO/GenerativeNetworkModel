function [energyout,ksout,netsout,ptsout] = ...
    fcn_sample_networks(A,Aseed,Dist,m,modeltype,modelvar,ndraw,nlvl,etalim,gamlim,pow)

indseed = find(triu(Aseed,1));
mseed = nnz(Aseed)/2;
mrem = m - mseed;
totalsamples = ndraw*nlvl;
ptsout = zeros(totalsamples,2);
energyout = zeros(totalsamples,1);
netsout = zeros(mrem,totalsamples);
ksout = zeros(totalsamples,4);

eta = unifrnd(etalim(1),etalim(2),ndraw,1);
gam = unifrnd(gamlim(1),gamlim(2),ndraw,1);
ptsout(1:ndraw,:) = [eta,gam];

fprintf('level %i of %i\n',1,nlvl);
[energyout(1:ndraw),netstemp,ksout(1:ndraw,:)] = fcn_eval_ks_max(modeltype,modelvar,m,Dist,A,Aseed,ptsout(1:ndraw,:));
for i = 1:mseed
    netstemp(netstemp == indseed(i)) = [];
end
netstemp = reshape(netstemp,[mrem,ndraw]);
netsout(:,1:ndraw) = netstemp;

powvals = linspace(0,pow,nlvl);
for ilvl = 2:nlvl
    
    pow = powvals(ilvl);
    fprintf('level %i of %i\n',ilvl,nlvl);
    ind = 1:ndraw*(ilvl - 1);
    
    if strcmp(modeltype,'sptl')
        ptsnew = fcn_voronoi_select_sptl(ptsout(ind,1),energyout(ind),ndraw,etalim,pow);
    else
        ptsnew = fcn_voronoi_select(ptsout(ind,:),energyout(ind),ndraw,etalim,gamlim,pow);
    end
    
    indnew = (1 + (ilvl - 1)*ndraw):(ilvl*ndraw);
    [energyout(indnew),netstemp,ksout(indnew,:)] = fcn_eval_ks_max(modeltype,modelvar,m,Dist,A,Aseed,ptsnew);
    if strcmp(modeltype,'sptl')
        ptsout(indnew,1) = ptsnew;
    else
        ptsout(indnew,:) = ptsnew;
    end
    for i = 1:mseed
        netstemp(netstemp == indseed(i)) = [];
    end
    netstemp = reshape(netstemp,[mrem,ndraw]);
    netsout(:,indnew) = netstemp;
    
end