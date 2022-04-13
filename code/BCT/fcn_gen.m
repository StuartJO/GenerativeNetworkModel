function b = fcn_gen(A,D,m,modeltype,modelvar,params)

% clear all
% close all
% clc
%
% load ../mats/griffa_mats_res3.mat
% A = A(:,:,1) >= 5;
% m = nnz(A)/2;
% D = mean(D,3);
% A = fcn_build_seed_network(A,D);
% modeltype = 'deg-min';
% params = [-3, 1; -1, 2; 0, 1; -3, 2];

n = length(D);
nparams = size(params,1);
mask = triu(ones(n),1) > 0;
[u,v] = find(mask);
indx = (v - 1)*n + u;
bseed = find(A(indx));
mseed = length(bseed);

switch modeltype
    
    case 'chain'
        b = zeros(m,nparams);
        for iparam = 1:nparams
            eta = params(iparam,1);
            gam = params(iparam,2);
            b(:,iparam) = fcn_chain(A,D,m,eta,gam,modelvar);
        end
        
    case 'chain2'
        b = zeros(m,nparams);
        for iparam = 1:nparams
            eta = params(iparam,1);
            gam = params(iparam,2);
            b(:,iparam) = fcn_chain2(A,D,m,eta,gam,modelvar);
        end
        
    case 'clu-avg'
        b = zeros(m,nparams);
        clu = clustering_coef_bu(A);
        Kseed = bsxfun(@plus,clu(:,ones(1,n)),clu')/2;
        for iparam = 1:nparams
            eta = params(iparam,1);
            gam = params(iparam,2);
            b(:,iparam) = fcn_clu_avg(A,Kseed,D,m,eta,gam,modelvar);
        end
        
    case 'clu2-avg'
        b = zeros(m,nparams);
        for iparam = 1:nparams
            eta = params(iparam,1);
            gam = params(iparam,2);
            b(:,iparam) = fcn_clu2avg(A,D,m,eta,gam);
        end
        
    case 'clu-diff'
        b = zeros(m,nparams);
        clu = clustering_coef_bu(A);
        Kseed = abs(bsxfun(@minus,clu(:,ones(1,n)),clu'));
        for iparam = 1:nparams
            eta = params(iparam,1);
            gam = params(iparam,2);
            b(:,iparam) = fcn_clu_max(A,Kseed,D,m,eta,gam,modelvar);
        end
        
    case 'clu2-diff'
        b = zeros(m,nparams);
        for iparam = 1:nparams
            eta = params(iparam,1);
            gam = params(iparam,2);
            b(:,iparam) = fcn_clu2diff(A,D,m,eta,gam);
        end
        
    case 'clu-max'
        b = zeros(m,nparams);
        clu = clustering_coef_bu(A);
        Kseed = bsxfun(@max,clu(:,ones(1,n)),clu');
        for iparam = 1:nparams
            eta = params(iparam,1);
            gam = params(iparam,2);
            b(:,iparam) = fcn_clu_max(A,Kseed,D,m,eta,gam,modelvar);
        end
        
    case 'clu2-max'
        b = zeros(m,nparams);
        for iparam = 1:nparams
            eta = params(iparam,1);
            gam = params(iparam,2);
            b(:,iparam) = fcn_clu2max(A,D,m,eta,gam);
        end
        
    case 'clu-min'
        b = zeros(m,nparams);
        clu = clustering_coef_bu(A);
        Kseed = bsxfun(@min,clu(:,ones(1,n)),clu');
        for iparam = 1:nparams
            eta = params(iparam,1);
            gam = params(iparam,2);
            b(:,iparam) = fcn_clu_min(A,Kseed,D,m,eta,gam,modelvar);
        end
        
    case 'clu2-min'
        b = zeros(m,nparams);
        for iparam = 1:nparams
            eta = params(iparam,1);
            gam = params(iparam,2);
            b(:,iparam) = fcn_clu2min(A,D,m,eta,gam);
        end
        
    case 'clu-prod'
        b = zeros(m,nparams);
        clu = clustering_coef_bu(A);
        Kseed = clu*clu';
        for iparam = 1:nparams
            eta = params(iparam,1);
            gam = params(iparam,2);
            b(:,iparam) = fcn_clu_prod(A,Kseed,D,m,eta,gam,modelvar);
        end
        
    case 'clu2-prod'
        b = zeros(m,nparams);
        for iparam = 1:nparams
            eta = params(iparam,1);
            gam = params(iparam,2);
            b(:,iparam) = fcn_clu2prod(A,D,m,eta,gam);
        end
        
    case 'deg-avg'
        b = zeros(m,nparams);
        kseed = sum(A,2);
        Kseed = bsxfun(@plus,kseed(:,ones(1,n)),kseed')/2;
        for iparam = 1:nparams
            eta = params(iparam,1);
            gam = params(iparam,2);
            b(:,iparam) = fcn_deg_avg(A,Kseed,D,m,eta,gam,modelvar);
        end
        
    case 'deg2-avg'
        b = zeros(m,nparams);
        for iparam = 1:nparams
            eta = params(iparam,1);
            gam = params(iparam,2);
            b(:,iparam) = fcn_deg2avg(A,D,m,eta,gam);
        end
        
    case 'deg-diff'
        b = zeros(m,nparams);
        kseed = sum(A,2);
        Kseed = abs(bsxfun(@minus,kseed(:,ones(1,n)),kseed'));
        for iparam = 1:nparams
            eta = params(iparam,1);
            gam = params(iparam,2);
            b(:,iparam) = fcn_deg_diff(A,Kseed,D,m,eta,gam,modelvar);
        end
        
    case 'deg2-diff'
        b = zeros(m,nparams);
        for iparam = 1:nparams
            eta = params(iparam,1);
            gam = params(iparam,2);
            b(:,iparam) = fcn_deg2diff(A,D,m,eta,gam);
        end
        
    case 'deg-max'
        b = zeros(m,nparams);
        kseed = sum(A,2);
        Kseed = bsxfun(@max,kseed(:,ones(1,n)),kseed');
        for iparam = 1:nparams
            eta = params(iparam,1);
            gam = params(iparam,2);
            b(:,iparam) = fcn_deg_max(A,Kseed,D,m,eta,gam,modelvar);
        end
        
    case 'deg2-max'
        b = zeros(m,nparams);
        for iparam = 1:nparams
            eta = params(iparam,1);
            gam = params(iparam,2);
            b(:,iparam) = fcn_deg2max(A,D,m,eta,gam);
        end
        
    case 'deg-min'
        b = zeros(m,nparams);
        b(1:mseed,:) = repmat(bseed,[1,nparams]);
        kseed = sum(A,2);
        Kseed = bsxfun(@min,kseed(:,ones(1,n)),kseed');
        for iparam = 1:nparams
            eta = params(iparam,1);
            gam = params(iparam,2);
            b(:,iparam) = fcn_deg_min(A,Kseed,D,m,eta,gam,modelvar);
        end
        
    case 'deg2-min'
        b = zeros(m,nparams);
        for iparam = 1:nparams
            eta = params(iparam,1);
            gam = params(iparam,2);
            b(:,iparam) = fcn_deg2min(A,D,m,eta,gam);
        end
        
    case 'deg-prod'
        b = zeros(m,nparams);
        kseed = sum(A,2);
        Kseed = (kseed*kseed').*~eye(n);
        for iparam = 1:nparams
            eta = params(iparam,1);
            gam = params(iparam,2);
            b(:,iparam) = fcn_deg_prod(A,Kseed,D,m,eta,gam,modelvar);
        end
        
    case 'deg2-prod'
        b = zeros(m,nparams);
        for iparam = 1:nparams
            eta = params(iparam,1);
            gam = params(iparam,2);
            b(:,iparam) = fcn_deg2prod(A,D,m,eta,gam);
        end
        
    case 'neighbors'
        b = zeros(m,nparams);
        Kseed = (A*A).*~eye(n);
        for iparam = 1:nparams
            eta = params(iparam,1);
            gam = params(iparam,2);
            b(:,iparam) = fcn_nghbrs(A,Kseed,D,m,eta,gam,modelvar);
        end
        
    case 'matching'
        b = zeros(m,nparams);
        Kseed = matching_ind(A);
        Kseed = Kseed + Kseed';
        for iparam = 1:nparams
            eta = params(iparam,1);
            gam = params(iparam,2);
            b(:,iparam) = fcn_matching(A,Kseed,D,m,eta,gam,modelvar);
        end
        
    case 'path-length'
        b = zeros(m,nparams);
        K = Soupor_FastFloyd(A);
        K = exp(K);
        for iparam = 1:nparams
            eta = params(iparam,1);
            gam = params(iparam,2);
            b(:,iparam) = fcn_path_length(A,K,D,m,eta,gam,modelvar);
        end
        
    case 'pref-sptl'
        b = zeros(m,nparams);
        for iparam = 1:nparams
            eta = params(iparam,1);
            gam = params(iparam,2);
            b(:,iparam) = fcn_pref_sptl(A,D,m,eta,gam,modelvar);
        end
        
    case 'sptl'
        if length(modelvar) > 1
            modelvar = modelvar{1};
        end
        b = zeros(m,nparams);
        for iparam = 1:nparams
            eta = params(iparam,1);
            b(:,iparam) = fcn_sptl(A,D,m,eta,modelvar);
        end
        
    case 'deg-avg-mod'
        Dnew = zeros(n);
        nci = max(ci);
        for u = 1:nci
            indu = ci == u;
            for v = 1:nci
                indv = ci == v;
                d = D(indu,indv);
                Dnew(indu,indv) = nanmean(d(:));
            end
        end
        D = Dnew;
        b = zeros(m,nparams);
        kseed = sum(A,2);
        Kseed = bsxfun(@plus,kseed(:,ones(1,n)),kseed')/2;
        for iparam = 1:nparams
            eta = params(iparam,1);
            gam = params(iparam,2);
            b(:,iparam) = fcn_deg_avg(A,Kseed,D,m,eta,gam,modelvar);
        end
        
    case 'deg-diff-mod'
        Dnew = zeros(n);
        nci = max(ci);
        for u = 1:nci
            indu = ci == u;
            for v = 1:nci
                indv = ci == v;
                d = D(indu,indv);
                Dnew(indu,indv) = nanmean(d(:));
            end
        end
        D = Dnew;
        b = zeros(m,nparams);
        kseed = sum(A,2);
        Kseed = abs(bsxfun(@minus,kseed(:,ones(1,n)),kseed'));
        for iparam = 1:nparams
            eta = params(iparam,1);
            gam = params(iparam,2);
            b(:,iparam) = fcn_deg_diff(A,Kseed,D,m,eta,gam,modelvar);
        end
        
    case 'deg-max-mod'
        Dnew = zeros(n);
        nci = max(ci);
        for u = 1:nci
            indu = ci == u;
            for v = 1:nci
                indv = ci == v;
                d = D(indu,indv);
                Dnew(indu,indv) = nanmean(d(:));
            end
        end
        D = Dnew;
        b = zeros(m,nparams);
        kseed = sum(A,2);
        Kseed = bsxfun(@max,kseed(:,ones(1,n)),kseed');
        for iparam = 1:nparams
            eta = params(iparam,1);
            gam = params(iparam,2);
            b(:,iparam) = fcn_deg_max(A,Kseed,D,m,eta,gam,modelvar);
        end
        
    case 'deg-min-mod'
        Dnew = zeros(n);
        nci = max(ci);
        for u = 1:nci
            indu = ci == u;
            for v = 1:nci
                indv = ci == v;
                d = D(indu,indv);
                Dnew(indu,indv) = nanmean(d(:));
            end
        end
        D = Dnew;
        b = zeros(m,nparams);
        b(1:mseed,:) = repmat(bseed,[1,nparams]);
        kseed = sum(A,2);
        Kseed = bsxfun(@min,kseed(:,ones(1,n)),kseed');
        for iparam = 1:nparams
            eta = params(iparam,1);
            gam = params(iparam,2);
            b(:,iparam) = fcn_deg_min(A,Kseed,D,m,eta,gam,modelvar);
        end
        
    case 'deg-prod-mod'
        Dnew = zeros(n);
        nci = max(ci);
        for u = 1:nci
            indu = ci == u;
            for v = 1:nci
                indv = ci == v;
                d = D(indu,indv);
                Dnew(indu,indv) = nanmean(d(:));
            end
        end
        D = Dnew;
        b = zeros(m,nparams);
        kseed = sum(A,2);
        Kseed = (kseed*kseed').*~eye(n);
        for iparam = 1:nparams
            eta = params(iparam,1);
            gam = params(iparam,2);
            b(:,iparam) = fcn_deg_prod(A,Kseed,D,m,eta,gam,modelvar);
        end
        
    case 'neighbors-mod'
        Dnew = zeros(n);
        nci = max(ci);
        for u = 1:nci
            indu = ci == u;
            for v = 1:nci
                indv = ci == v;
                d = D(indu,indv);
                Dnew(indu,indv) = nanmean(d(:));
            end
        end
        D = Dnew;
        b = zeros(m,nparams);
        Kseed = (A*A).*~eye(n);
        for iparam = 1:nparams
            eta = params(iparam,1);
            gam = params(iparam,2);
            b(:,iparam) = fcn_nghbrs(A,Kseed,D,m,eta,gam,modelvar);
        end
        
    case 'pref-sptl-mod'
        Dnew = zeros(n);
        nci = max(ci);
        for u = 1:nci
            indu = ci == u;
            for v = 1:nci
                indv = ci == v;
                d = D(indu,indv);
                Dnew(indu,indv) = nanmean(d(:));
            end
        end
        D = Dnew;
        b = zeros(m,nparams);
        for iparam = 1:nparams
            eta = params(iparam,1);
            gam = params(iparam,2);
            b(:,iparam) = fcn_pref_sptl(A,D,m,eta,gam,modelvar);
        end
        
    case 'sptl-mod'
        Dnew = zeros(n);
        nci = max(ci);
        for u = 1:nci
            indu = ci == u;
            for v = 1:nci
                indv = ci == v;
                d = D(indu,indv);
                Dnew(indu,indv) = nanmean(d(:));
            end
        end
        D = Dnew;
        if length(modelvar) > 1
            modelvar = modelvar{1};
        end
        b = zeros(m,nparams);
        for iparam = 1:nparams
            eta = params(iparam,1);
            b(:,iparam) = fcn_sptl(A,D,m,eta,modelvar);
        end
        
end
