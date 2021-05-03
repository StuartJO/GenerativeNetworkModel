function [B,b] = gen_model_additive_fast2_PD_normsum(A,PDMs,m,modeltype,modelvar,eta,gam,alpha,epsilon)
% GENERATIVE_MODEL          Run generative model code
%
%   B = GENERATIVE_MODEL(A,D,m,modeltype,modelvar,params)
%
%   Generates synthetic networks using the models described in the study by
%   Betzel et al (2016) in Neuroimage.
%
%   Inputs:
%           A,          binary network of seed connections
%           D,          Euclidean distance/fiber length matrix
%           O,          Other distance metrics, with each as the third
%                       dimension in a matrix
%           m,          number of connections that should be present in
%                       final synthetic network
%           modeltype,  specifies the generative rule (see below)
%           modelvar,   specifies whether the generative rules are based on
%                       power-law or exponential relationship
%                       ({'powerlaw'}|{'exponential})
%           params,     either a vector (in the case of the geometric
%                       model) or a matrix (for all other models) of
%                       parameters at which the model should be evaluated.
%           Oparams,    a vector of parameters to use for the other
%                       distance metrics where Oparams(i) is the parameter
%                       to use for O{i}
%           epsilon,    the baseline probability of forming a particular
%                       connection (should be a very small number
%                       {default = 1e-5}).
%
%   Output:
%           B,          an adjacency network
%
%
%   Full list of model types:
%   (each model type realizes a different generative rule)
%
%       1.  'sptl'          spatial model
%       2.  'neighbors'     number of common neighbors
%       3.  'matching'      matching index
%       4.  'clu-avg'       average clustering coeff.
%       5.  'clu-min'       minimum clustering coeff.
%       6.  'clu-max'       maximum clustering coeff.
%       7.  'clu-diff'      difference in clustering coeff.
%       8.  'clu-prod'      product of clustering coeff.
%       9.  'deg-avg'       average degree
%       10. 'deg-min'       minimum degree
%       11. 'deg-max'       maximum degree
%       12. 'deg-diff'      difference in degree
%       13. 'deg-prod'      product of degree
%       14. 'com'           communicability
%
%   Example usage:
%
%       load demo_generative_models_data
%
%       % get number of bi-directional connections
%       m = nnz(A)/2;
%
%       % get cardinality of network
%       n = length(A);
%
%       % set model type
%       modeltype = 'neighbors';
%
%       % set whether the model is based on powerlaw or exponentials
%       modelvar = [{'powerlaw'},{'powerlaw'}];
%
%       % choose some model parameters
%       params = [-2,0.2; -5,1.2; -1,1.5];
%       nparams = size(params,1);
%
%       % generate synthetic networks
%       B = generative_model(Aseed,D,m,modeltype,modelvar,params);
%
%       % store them in adjacency matrix format
%       Asynth = zeros(n,n,nparams);
%       for i = 1:nparams;
%           a = zeros(n); a(B(:,i)) = 1; a = a + a';
%           Asynth(:,:,i) = a;
%       end
%
%   Reference: Betzel et al (2016) Neuroimage 124:1054-64.
%
%   Richard Betzel, Indiana University/University of Pennsylvania, 2015

if ~exist('epsilon','var')
    epsilon = 0;
end



if iscell(PDMs) || size(PDMs,3) > 1
    
    if iscell(PDMs)
        if length(PDMs) == 2
             PD = PDMs{2};
             D = PDMs{1};
        else
            PD = [];
            D = PDMs{1};
        end
    else
          PD = squeeze(PDMs(:,:,2));
            D = squeeze(PDMs(:,:,1));
          
    end
else
            PD = [];
            D = PDMs;

end
 n = length(D);
if isempty(PD)
    
    PD = ones(n);
    alpha(3) = 0;
    eta(2) = 1;   
end

% if ~exist('PD','var') && size(D,3) == 1 && ~iscell(D)
%     
%     PD = ones(n);
%     alpha(3) = 0;
%     eta(2) = 1;
% 
% elseif ~exist('PD','var') && iscell(D)
%     PD = D{2};
%     D = D{1};
% elseif ~exist('PD','var') && size(D,3) == 1
%     PD = squeeze(D(:,:,2));
%     D = squeeze(D(:,:,1));
% end

switch modeltype

    case 'clu-avg'
        clu = clustering_coef_bu(A);
        Kseed = bsxfun(@plus,clu(:,ones(1,n)),clu')/2;
        b = fcn_clu(modeltype,A,Kseed,D,m,eta,gam,modelvar,epsilon,PD,alpha);

    case 'clu-diff'
        clu = clustering_coef_bu(A);
        Kseed = abs(bsxfun(@minus,clu(:,ones(1,n)),clu'));
        b = fcn_clu(modeltype,A,Kseed,D,m,eta,gam,modelvar,epsilon,PD,alpha);

    case 'clu-max'
        clu = clustering_coef_bu(A);
        Kseed = bsxfun(@max,clu(:,ones(1,n)),clu');
        b = fcn_clu(modeltype,A,Kseed,D,m,eta,gam,modelvar,epsilon,PD,alpha);

    case 'clu-min'
        clu = clustering_coef_bu(A);
        Kseed = bsxfun(@min,clu(:,ones(1,n)),clu');
        b = fcn_clu(modeltype,A,Kseed,D,m,eta,gam,modelvar,epsilon,PD,alpha);

    case 'clu-prod'
        clu = clustering_coef_bu(A);
        Kseed = clu*clu';
        b = fcn_clu(modeltype,A,Kseed,D,m,eta,gam,modelvar,epsilon,PD,alpha);

    case 'deg-avg'
        kseed = sum(A,2);
        Kseed = bsxfun(@plus,kseed(:,ones(1,n)),kseed')/2;
        b = fcn_deg(modeltype,A,Kseed,D,m,eta,gam,modelvar,epsilon,PD,alpha);

    case 'deg-diff'
        kseed = sum(A,2);
        Kseed = abs(bsxfun(@minus,kseed(:,ones(1,n)),kseed'));
        b = fcn_deg(modeltype,A,Kseed,D,m,eta,gam,modelvar,epsilon,PD,alpha);

    case 'deg-max'
        kseed = sum(A,2);
        Kseed = bsxfun(@max,kseed(:,ones(1,n)),kseed');
        b = fcn_deg(modeltype,A,Kseed,D,m,eta,gam,modelvar,epsilon,PD,alpha);

    case 'deg-min'
        kseed = sum(A,2);
        Kseed = bsxfun(@min,kseed(:,ones(1,n)),kseed');
        b = fcn_deg(modeltype,A,Kseed,D,m,eta,gam,modelvar,epsilon,PD,alpha);

    case 'deg-prod'
        kseed = sum(A,2);
        Kseed = (kseed*kseed').*~eye(n);
        b = fcn_deg(modeltype,A,Kseed,D,m,eta,gam,modelvar,epsilon,PD,alpha);

    case 'neighbors'
        Kseed = (A*A).*~eye(n);
        b = fcn_nghbrs(A,Kseed,D,m,eta,gam,modelvar,epsilon,PD,alpha);

    case 'matching'
        Kseed = matching_ind(A);
        Kseed = Kseed + Kseed';
        b = fcn_matching(A,Kseed,D,m,eta,gam,modelvar,epsilon,PD,alpha);

    case 'sptl'
        b = fcn_sptl(A,D,m,eta,modelvar,PD,alpha);

    case 'com'
        Kseed = gexpm(A);
        b = fcn_com(A,Kseed,D,m,eta,gam,modelvar,epsilon,PD,alpha);
end

B = zeros(n);
B(b) = 1;
B = B + B';

function b = fcn_clu(modeltype,A,K,D,m,eta,gam,modelvar,epsilon,PD,alpha)
K = K + epsilon;

n = length(D);
mseed = nnz(A)/2;
A = A > 0;
mv1 = modelvar{1};
mv2 = modelvar{2};
mv3 = modelvar{3};

[u,v] = find(triu(ones(n),1));
indx = (v - 1)*n + u;

EdgesTriu = zeros(size(indx));
b = zeros(m,1);
b(1:mseed) = find(A(indx));

EdgesTriu(b(b~=0)) = 1;

switch mv1
    case 'powerlaw'
        Df = D.^eta(1);
    case 'exponential'       
        Df = exp(eta(1)*(D));    
end
switch mv3
    case 'powerlaw'
        PDf = PD.^eta(2);
    case 'exponential'       
        PDf = exp(eta(2)*(PD));    
end

Fd = alpha(1)*(Df./sum(Df(indx(~EdgesTriu)))) + alpha(3)*(PDf./sum(PDf(indx(~EdgesTriu))));

switch mv2
    case 'powerlaw'
        Kf = K.^gam;
    case 'exponential'       
        Kf = exp(gam*K);    
end
Kf(isinf(Kf)) = 0;
Fk = alpha(2)*(Kf./sum(Kf(indx(~EdgesTriu))));
Fk(isnan(Fk)) = 0;Fk(isinf(Fk)) = 0;

c = clustering_coef_bu(A);
k = sum(A,2);

Ff = (Fd+Fk).*~A;

P = Ff(indx);

for i = (mseed + 1):m
    C = [0; cumsum(P)];
    r = sum(rand*C(end) >= C);
    
    b(i) = r;
    EdgesTriu(b(i)) = 1;
    uu = u(r);
    vv = v(r);
    A(uu,vv) = 1;
    A(vv,uu) = 1;
    k([uu,vv]) = k([uu,vv]) + 1;
    bu = A(uu,:);
    su = A(bu,bu);
    bv = A(vv,:);
    sv = A(bv,bv);
    bth = bu & bv;
    c(bth) = c(bth) + 2./(k(bth).^2 - k(bth));
    c(uu) = nnz(su)/(k(uu)*(k(uu) - 1));
    c(vv) = nnz(sv)/(k(vv)*(k(vv) - 1));
    c(k <= 1) = 0;
    bth([uu,vv]) = true;
    
    switch modeltype
        case 'clu-avg'
    
            K(:,bth) = bsxfun(@plus,c(:,ones(1,sum(bth))),c(bth,:)')/2 + epsilon;
            K(bth,:) = bsxfun(@plus,c(:,ones(1,sum(bth))),c(bth,:)')'/2 + epsilon;

        case 'clu-diff'

            K(:,bth) = abs(bsxfun(@minus,c(:,ones(1,sum(bth))),c(bth,:)')) + epsilon;
            K(bth,:) = abs(bsxfun(@minus,c(:,ones(1,sum(bth))),c(bth,:)'))' + epsilon;

        case 'clu-max'
            
            K(:,bth) = bsxfun(@max,c(:,ones(1,sum(bth))),c(bth,:)') + epsilon;
            K(bth,:) = bsxfun(@max,c(:,ones(1,sum(bth))),c(bth,:)')' + epsilon;
            
        case 'clu-min'
            
            K(:,bth) = bsxfun(@min,c(:,ones(1,sum(bth))),c(bth,:)') + epsilon;
            K(bth,:) = bsxfun(@min,c(:,ones(1,sum(bth))),c(bth,:)')' + epsilon;
    
        case 'clu-prod'
            
            K(bth,:) = (c(bth,:)*c') + epsilon;
            K(:,bth) = (c*c(bth,:)') + epsilon;
    end
    
    Fd = alpha(1)*(Df./sum(Df(indx(~EdgesTriu)))) + alpha(3)*(PDf./sum(PDf(indx(~EdgesTriu))));
    
    switch mv2
        case 'powerlaw'
            Kf = K.^gam;
        case 'exponential'       
            Kf = exp(gam*K);    
    end
    Kf(isinf(Kf)) = 0;
    Fk = alpha(2)*(Kf./sum(sum(Kf.*~A)));
    Fk(isnan(Fk)) = 0;Fk(isinf(Fk)) = 0;
    Ff = (Fd+Fk).*~A;
    P = Ff(indx);
    
end
b = indx(b);


function b = fcn_deg(modeltype,A,K,D,m,eta,gam,modelvar,epsilon,PD,alpha)
n = length(D);
mseed = nnz(A)/2;
k = sum(A,2);
[u,v] = find(triu(ones(n),1));
indx = (v - 1)*n + u;

EdgesTriu = zeros(size(indx));
b = zeros(m,1);
b(1:mseed) = find(A(indx));

EdgesTriu(b(b~=0)) = 1;
D = D(indx);
PD = PD(indx);
mv1 = modelvar{1};
mv2 = modelvar{2};
mv3 = modelvar{3};

switch mv1
    case 'powerlaw'
        Df = D.^eta(1);
    case 'exponential'       
        Df = exp(eta(1)*(D));    
end
switch mv3
    case 'powerlaw'
        PDf = PD.^eta(2);
    case 'exponential'       
        PDf = exp(eta(2)*(PD));    
end

Fd = alpha(1)*(Df./sum(Df(~EdgesTriu))) + alpha(3)*(PDf./sum(PDf(~EdgesTriu)));
K = K + epsilon;

switch mv2
    case 'powerlaw'
        Kf = K.^gam;
    case 'exponential'       
        Kf = exp(gam*K);    
end
Kf(isinf(Kf)) = 0;
Fk = alpha(2)*(Kf./sum(Kf(indx(~EdgesTriu))));
Fk(isnan(Fk)) = 0;Fk(isinf(Fk)) = 0;
    
    
P = (Fd+Fk(indx)).*~A(indx);

for i = (mseed + 1):m
    C = [0; cumsum(P)];
    r = sum(rand*C(end) >= C);
    w = [u(r),v(r)];
    k(w) = k(w) + 1;
    
    
    switch modeltype
        case 'deg-avg'
    
            K(:,w) = [((k + k(w(1)))/2) + epsilon, ((k + k(w(2)))/2) + epsilon];
            K(w,:) = ([((k + k(w(1)))/2) + epsilon, ((k + k(w(2)))/2) + epsilon])';
    
        case 'deg-diff'
    
            K(:,w) = (abs([k - k(w(1)), k - k(w(2))]) + epsilon);
            K(w,:) = (abs([k - k(w(1)), k - k(w(2))]) + epsilon)';
    
        case 'deg-min'
            
            K(:,w) = [min(k,k(w(1))) + epsilon, min(k,k(w(2))) + epsilon];
            K(w,:) = [min(k,k(w(1))) + epsilon, min(k,k(w(2))) + epsilon]'; 
            
            
        case 'deg-max'
            
            K(:,w) = [max(k,k(w(1))) + epsilon, max(k,k(w(2))) + epsilon];
            K(w,:) = [max(k,k(w(1))) + epsilon, max(k,k(w(2))) + epsilon]'; 
            
        case 'deg-prod'
            
            K(:,w) = [k*k(w(1)) + epsilon, k*k(w(2)) + epsilon];
            K(w,:) = [k*k(w(1)) + epsilon, k*k(w(2)) + epsilon]';
    
    end
    
    b(i) = r;
    EdgesTriu(b(i)) = 1;
    
    Fd = alpha(1)*(Df./sum(Df(~EdgesTriu))) + alpha(3)*(PDf./sum(PDf(~EdgesTriu)));
    
    switch mv2
        case 'powerlaw'
            Kf = K.^gam;
        case 'exponential'       
            Kf = exp(gam*K);    
    end
    Kf(isinf(Kf)) = 0;
    Fk = alpha(2)*(Kf./sum(Kf(indx(~EdgesTriu))));
    Fk(isnan(Fk)) = 0;Fk(isinf(Fk)) = 0;
    P = Fd+Fk(indx);

    P(b(1:i)) = 0;
    
end
b = indx(b);

function b = fcn_nghbrs(A,K,D,m,eta,gam,modelvar,epsilon,PD,alpha)
K = K + epsilon;

n = length(D);
mseed = nnz(A)/2;
A = A > 0;
mv1 = modelvar{1};
mv2 = modelvar{2};
mv3 = modelvar{3};

[u,v] = find(triu(ones(n),1));
indx = (v - 1)*n + u;

EdgesTriu = zeros(size(indx));
b = zeros(m,1);
b(1:mseed) = find(A(indx));

EdgesTriu(b(b~=0)) = 1;
switch mv1
    case 'powerlaw'
        Df = D.^eta(1);
    case 'exponential'       
        Df = exp(eta(1)*(D));    
end
switch mv3
    case 'powerlaw'
        PDf = PD.^eta(2);
    case 'exponential'       
        PDf = exp(eta(2)*(PD));    
end
Fd = alpha(1)*(Df./sum(Df(indx(~EdgesTriu)))) + alpha(3)*(PDf./sum(PDf(indx(~EdgesTriu))));
switch mv2
    case 'powerlaw'
        Kf = K.^gam;
    case 'exponential'       
        Kf = exp(gam*K);    
end
Kf(isinf(Kf)) = 0;
Fk = alpha(2)*(Kf./sum(Kf(indx(~EdgesTriu))));
Fk(isnan(Fk)) = 0;Fk(isinf(Fk)) = 0;
Ff = (Fd+Fk).*~A;

P = Ff(indx);


for i = (mseed + 1):m
    C = [0; cumsum(P)];
    r = sum(rand*C(end) >= C);
    b(i) = r;
    EdgesTriu(b(i)) = 1;
    uu = u(r);
    vv = v(r);
    x = A(uu,:);
    y = A(:,vv);
    A(uu,vv) = 1;
    A(vv,uu) = 1;
    K(uu,y) = K(uu,y) + 1;
    K(y,uu) = K(y,uu) + 1;
    K(vv,x) = K(vv,x) + 1;
    K(x,vv) = K(x,vv) + 1;
    
    Fd = alpha(1)*(Df./sum(Df(indx(~EdgesTriu)))) + alpha(3)*(PDf./sum(PDf(indx(~EdgesTriu))));
    switch mv2
        case 'powerlaw'
            Kf = K.^gam;
        case 'exponential'       
            Kf = exp(gam*K);    
    end
    Kf(isinf(Kf)) = 0;
    Fk = alpha(2)*(Kf./sum(sum(Kf.*~A)));
    Fk(isnan(Fk)) = 0;Fk(isinf(Fk)) = 0;
    
    Ff = Fd + Fk;

    Ff(A) = 0;

    P = Ff(indx);
    
end
b = indx(b);

function b = fcn_matching(A,K,D,m,eta,gam,modelvar,epsilon,PD,alpha)
K = K + epsilon;

n = length(D);
mseed = nnz(A)/2;
mv1 = modelvar{1};
mv2 = modelvar{2};
mv3 = modelvar{3};

[u,v] = find(triu(ones(n),1));
indx = (v - 1)*n + u;

EdgesTriu = zeros(size(indx));
b = zeros(m,1);
b(1:mseed) = find(A(indx));

EdgesTriu(b(b~=0)) = 1;
switch mv1
    case 'powerlaw'
        Df = D.^eta(1);
    case 'exponential'       
        Df = exp(eta(1)*(D));    
end
switch mv3
    case 'powerlaw'
        PDf = PD.^eta(2);
    case 'exponential'       
        PDf = exp(eta(2)*(PD));    
end
Fd = alpha(1)*(Df./sum(Df(indx(~EdgesTriu)))) + alpha(3)*(PDf./sum(PDf(indx(~EdgesTriu))));
switch mv2
    case 'powerlaw'
        Kf = K.^gam;
    case 'exponential'       
        Kf = exp(gam*K);    
end
Kf(isinf(Kf)) = 0;
Fk = alpha(2)*(Kf./sum(Kf(indx(~EdgesTriu))));
Fk(isnan(Fk)) = 0;Fk(isinf(Fk)) = 0;
Ff = (Fd+Fk).*~A;

P = Ff(indx);

for ii = (mseed + 1):m

    C = [0; cumsum(P)];
    r = sum(rand*C(end) >= C);
    
    b(ii) = r;
%     if r == 0
%        save('ERRORTEST.mat','C','r','P','Ff','Fk','EdgesTriu','ii','Fd','Kf') 
%     end
    EdgesTriu(b(ii)) = 1;
    uu = u(r);
    vv = v(r);

    A(uu,vv) = 1;
    A(vv,uu) = 1;

    updateuu = find(A*A(:,uu));
    updateuu(updateuu == uu) = [];
    updateuu(updateuu == vv) = [];

    updatevv = find(A*A(:,vv));
    updatevv(updatevv == uu) = [];
    updatevv(updatevv == vv) = [];

    c1 = [A(:,uu)', A(uu,:)];
    for i = 1:length(updateuu)
        j = updateuu(i);
        c2 = [A(:,j)' A(j,:)];
        use = ~(~c1&~c2);
        use(uu) = 0;  use(uu+n) = 0;
        use(j) = 0;  use(j+n) = 0;
        ncon = sum(c1(use))+sum(c2(use));
        if (ncon==0)
            K(uu,j) = epsilon;
            K(j,uu) = epsilon;
        else
            K(uu,j) = (2*(sum(c1(use)&c2(use))/ncon)) + epsilon;
            K(j,uu) = K(uu,j);
        end

    end

    c1 = [A(:,vv)', A(vv,:)];
    for i = 1:length(updatevv)
        j = updatevv(i);
        c2 = [A(:,j)' A(j,:)];
        use = ~(~c1&~c2);
        use(vv) = 0;  use(vv+n) = 0;
        use(j) = 0;  use(j+n) = 0;
        ncon = sum(c1(use))+sum(c2(use));
        if (ncon==0)
            K(vv,j) = epsilon;
            K(j,vv) = epsilon;
        else
            K(vv,j) = (2*(sum(c1(use)&c2(use))/ncon)) + epsilon;
            K(j,vv) = K(vv,j);
        end
    end
    Fd = alpha(1)*(Df./sum(Df(indx(~EdgesTriu)))) + alpha(3)*(PDf./sum(PDf(indx(~EdgesTriu))));
    switch mv2
        case 'powerlaw'
            Kf = K.^gam;
        case 'exponential'       
            Kf = exp(gam*K);    
    end
    Kf(isinf(Kf)) = 0;
    Fk = alpha(2)*(Kf./sum(sum(Kf.*~A)));
    Fk(isnan(Fk)) = 0;Fk(isinf(Fk)) = 0;
    Ff = (Fd+Fk).*~A;
    P = Ff(indx);
end
b = indx(b);

function b = fcn_sptl(A,D,m,eta,modelvar,PD,alpha)
n = length(D);
mseed = nnz(A)/2;
[u,v] = find(triu(ones(n),1));
indx = (v - 1)*n + u;

mv1 = modelvar{1};
mv3 = modelvar{3};

EdgesTriu = zeros(size(indx));
b = zeros(m,1);
b(1:mseed) = find(A(indx));

EdgesTriu(b(b~=0)) = 1;
switch mv1
    case 'powerlaw'
        Df = D.^eta(1);
    case 'exponential'       
        Df = exp(eta(1)*(D));    
end
switch mv3
    case 'powerlaw'
        PDf = PD.^eta(2);
    case 'exponential'       
        PDf = exp(eta(2)*(PD));    
end
Fd = alpha(1)*(Df./sum(Df(indx(~EdgesTriu)))) + alpha(3)*(PDf./sum(PDf(indx(~EdgesTriu))));

P = Fd(indx).*~A(indx);

for i = (mseed + 1):m
    C = [0; cumsum(P)];
    r = sum(rand*C(end) >= C);
    b(i) = r;

    EdgesTriu(b(i)) = 1;
    
    Fd = alpha(1)*(Df./sum(Df(indx(~EdgesTriu)))) + alpha(3)*(PDf./sum(PDf(indx(~EdgesTriu))));
    
    P = Fd(indx);
    P(b(1:i)) = 0;
end
b = indx(b);

function b = fcn_com(A,K,D,m,eta,gam,modelvar,epsilon,PD,alpha)

K = K + epsilon;

n = length(D);
mseed = nnz(A)/2;
mv1 = modelvar{1};
mv2 = modelvar{2};
mv3 = modelvar{3};

[u,v] = find(triu(ones(n),1));
indx = (v - 1)*n + u;

EdgesTriu = zeros(size(indx));
b = zeros(m,1);
b(1:mseed) = find(A(indx));

EdgesTriu(b(b~=0)) = 1;
switch mv1
    case 'powerlaw'
        Df = D.^eta(1);
    case 'exponential'       
        Df = exp(eta(1)*(D));    
end
switch mv3
    case 'powerlaw'
        PDf = PD.^eta(2);
    case 'exponential'       
        PDf = exp(eta(2)*(PD));    
end
Fd = alpha(1)*(Df./sum(Df(indx(~EdgesTriu)))) + + alpha(3)*(PDf./sum(PDf(indx(~EdgesTriu))));
switch mv2
    case 'powerlaw'
        Kf = K.^gam;
    case 'exponential'       
        Kf = exp(gam*K);    
end
Kf(isinf(Kf)) = 0;
Fk = alpha(2)*(Kf./sum(Kf(indx(~EdgesTriu))));
Fk(isnan(Fk)) = 0;Fk(isinf(Fk)) = 0;
Ff = (Fd+Fk).*~A;

P = Ff(indx);
%size(P)

for i = (mseed + 1):m
    C = [0; cumsum(P)];
    r = sum(rand*C(end) >= C);
    
    b(i) = r;
    EdgesTriu(b(i)) = 1;
    %r
    uu = u(r);
    vv = v(r);
    A(uu,vv) = 1;
    A(vv,uu) = 1;
    K = gexpm(A)+epsilon;
    Fd = alpha(1)*(Df./sum(Df(indx(~EdgesTriu)))) + alpha(3)*(PDf./sum(PDf(indx(~EdgesTriu))));
    switch mv2
        case 'powerlaw'
            Kf = K.^gam;
        case 'exponential'       
            Kf = exp(gam*(K));    
    end
    Kf(isinf(Kf)) = 0;
    Fk = alpha(2)*(Kf./sum(sum(Kf.*~A)));
    Fk(isnan(Fk)) = 0;Fk(isinf(Fk)) = 0;
    Ff = (Fd+Fk).*~A;
    P = Ff(indx);
end
b = indx(b);