function [B,b] = gen_model_mult(A,PD,m,modeltype,modelvar,PDexpo,gam,epsilon)
% gen_model_mult          Run generative model code for the multiplicative
% model
%
%   Generates synthetic networks using the models described in the study by
%   Betzel et al (2016) in Neuroimage.
%
%   Inputs:
%           A,          binary network of seed connections
%           PD,         Euclidean distance/fiber length/node similarity
%                       matrix. Multiple can be input either as a cell, 
%                       where each cell contains a different matrix or as a
%                       3D matrix (n*n*nPD, where n is the number of nodes
%                       and nPD is the number of PD matrices).
%           m,          number of connections that should be present in
%                       final synthetic network
%           modeltype,  specifies the generative rule (see below)
%           modelvar,   specifies whether the generative rules are based on
%                       power-law or exponential relationship
%                       ({'powerlaw'}|{'exponential})
%           PDexpo,     the parameter controlling the values in PD. If
%                       there are multipe PD matrices, PDexpo should be a
%                       vector where each index gives the marameter for the
%                       corresponding PD matrix
%           gam,        the parameter controlling topology
%           epsilon,    the baseline probability of forming a particular
%                       connection (should be a very small number
%                       {default = 1e-6}).
%
%   Output:
%           B,          an adjacency matrix
%           b,          a vector giving the index of each edge in B. Note
%                       that the ordering of b shows which edges formed
%                       first (e.g., b(1) was the fiorst edge to form, b(2)
%                       the second etc etc).
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
%       How to convert b to B:
%       n = length(A); B = zeros(n); B(b(:,i)) = 1; B = B + B'; 
%
%
%   Reference: Betzel et al (2016) Neuroimage 124:1054-64.
%
%   Richard Betzel, Indiana University/University of Pennsylvania, 2015
%   Edited by Stuart Oldham, Monash University 2021

if ~exist('epsilon','var')
    epsilon = 1e-6;
end

n = length(A);

% Perform the multiplication of PDs as these values will not change across
% iterations

nPD = length(PD);

Df = zeros(n,n,nPD);

mv1 = modelvar{1};

if iscell(mv1)
    for ii = 1:nPD
         switch mv1{ii}
            case 'powerlaw'
                    Df(:,:,ii) = PD{ii}.^PDexpo(ii);
            case 'exponential'       
                    Df(:,:,ii) = exp(PDexpo(ii)*(PD{ii}));
         end    
    end    
else
    switch mv1
        case 'powerlaw'
            for i = 1:nPD
                Df(:,:,i) = PD{i}.^PDexpo(i);
            end
        case 'exponential'       
            for i = 1:nPD
                Df(:,:,i) = exp(PDexpo(i)*(PD{i}));  
            end
     end
end

Fd = prod(Df,3);

switch modeltype
    
    case 'clu-avg'
        clu = clustering_coef_bu(A);
        Kseed = bsxfun(@plus,clu(:,ones(1,n)),clu')/2;
        b = fcn_clu(modeltype,A,Kseed,Fd,m,gam,modelvar,epsilon);
        
    case 'clu-diff'
        clu = clustering_coef_bu(A);
        Kseed = abs(bsxfun(@minus,clu(:,ones(1,n)),clu'));
        b = fcn_clu(modeltype,A,Kseed,Fd,m,gam,modelvar,epsilon);
        
    case 'clu-max'
        clu = clustering_coef_bu(A);
        Kseed = bsxfun(@max,clu(:,ones(1,n)),clu');
        b = fcn_clu(modeltype,A,Kseed,Fd,m,gam,modelvar,epsilon);
        
    case 'clu-min'
        clu = clustering_coef_bu(A);
        Kseed = bsxfun(@min,clu(:,ones(1,n)),clu');
        b = fcn_clu(modeltype,A,Kseed,Fd,m,gam,modelvar,epsilon);
        
    case 'clu-prod'
        clu = clustering_coef_bu(A);
        Kseed = clu*clu';
        b = fcn_clu(modeltype,A,Kseed,Fd,m,gam,modelvar,epsilon);
        
    case 'deg-avg'
        kseed = sum(A,2);
        Kseed = bsxfun(@plus,kseed(:,ones(1,n)),kseed')/2;
        b = fcn_deg(modeltype,A,Kseed,Fd,m,gam,modelvar,epsilon);
        
    case 'deg-diff'
        kseed = sum(A,2);
        Kseed = abs(bsxfun(@minus,kseed(:,ones(1,n)),kseed'));
        b = fcn_deg(modeltype,A,Kseed,Fd,m,gam,modelvar,epsilon);
        
    case 'deg-max'
        kseed = sum(A,2);
        Kseed = bsxfun(@max,kseed(:,ones(1,n)),kseed');
        b = fcn_deg(modeltype,A,Kseed,Fd,m,gam,modelvar,epsilon);
        
    case 'deg-min'
        kseed = sum(A,2);
        Kseed = bsxfun(@min,kseed(:,ones(1,n)),kseed');
        b = fcn_deg(modeltype,A,Kseed,Fd,m,gam,modelvar,epsilon);
        
    case 'deg-prod'
        kseed = sum(A,2);
        Kseed = (kseed*kseed').*~eye(n);
        b = fcn_deg(modeltype,A,Kseed,Fd,m,gam,modelvar,epsilon);
        
    case 'neighbors'
        Kseed = (A*A).*~eye(n);
        b = fcn_nghbrs(A,Kseed,Fd,m,gam,modelvar,epsilon);
        
    case 'matching'
        Kseed = matching_ind(A);
        Kseed = Kseed + Kseed';
        b = fcn_matching(A,Kseed,Fd,m,gam,modelvar,epsilon);
        
    case 'sptl'
        b = fcn_sptl(A,Fd,m);
        
    case 'com'
        Kseed = gexpm(A);
        b = fcn_com(A,Kseed,Fd,m,gam,modelvar,epsilon);
end

B = zeros(n);
B(b) = 1;
B = B + B';

function b = fcn_clu(modeltype,A,K,Fd,m,gam,modelvar,epsilon)
K = K + epsilon;
n = length(Fd);
mseed = nnz(A)/2;
A = A > 0;

mv2 = modelvar{2};

switch mv2
    case 'powerlaw'
        Fk = K.^gam;
    case 'exponential'
        Fk = exp(gam*K);
end

c = clustering_coef_bu(A);
k = sum(A,2);

Ff = Fd.*Fk.*~A;
[u,v] = find(triu(ones(n),1));
indx = (v - 1)*n + u;
P = Ff(indx);
b = zeros(m,1);
b(1:mseed) = find(A(indx));
for i = (mseed + 1):m
    C = [0; cumsum(P)];
    r = sum(rand*C(end) >= C);
    b(i) = r;
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

    switch mv2
        case 'powerlaw'
            Ff(bth,:) = Fd(bth,:).*((K(bth,:)).^gam);
            Ff(:,bth) = Fd(:,bth).*((K(:,bth)).^gam);
        case 'exponential'
            Ff(bth,:) = Fd(bth,:).*exp((K(bth,:))*gam);
            Ff(:,bth) = Fd(:,bth).*exp((K(:,bth))*gam);
    end
    Ff = Ff.*~A;
    P = Ff(indx);
end
b = indx(b);

function b = fcn_deg(modeltype,A,K,Fd,m,gam,modelvar,epsilon)
n = length(Fd);
mseed = nnz(A)/2;
k = sum(A,2);
[u,v] = find(triu(ones(n),1));
indx = (v - 1)*n + u;
Fd = Fd(indx);

mv2 = modelvar{2};

K = K + epsilon;
switch mv2
    case 'powerlaw'
        Fk = K.^gam;
    case 'exponential'
        Fk = exp(gam*K);
end
P = Fd.*Fk(indx).*~A(indx);
b = zeros(m,1);
b(1:mseed) = find(A(indx));
for i = (mseed + 1):m
    C = [0; cumsum(P)];
    r = sum(rand*C(end) >= C);
    b(i) = r;
    w = [u(r),v(r)];
    k(w) = k(w) + 1;
    
    
    switch modeltype
        
    case 'deg-avg'
        
        switch mv2
        case 'powerlaw'
            Fk(:,w) = [((k + k(w(1)))/2) + epsilon, ((k + k(w(2)))/2) + epsilon].^gam;
            Fk(w,:) = ([((k + k(w(1)))/2) + epsilon, ((k + k(w(2)))/2) + epsilon].^gam)';
        case 'exponential'
            Fk(:,w) = exp([((k + k(w(1)))/2) + epsilon, ((k + k(w(2)))/2) + epsilon]*gam);
            Fk(w,:) = exp([((k + k(w(1)))/2) + epsilon, ((k + k(w(2)))/2) + epsilon]*gam)';
        end    

    case 'deg-diff'
        switch mv2
        case 'powerlaw'
            Fk(:,w) = (abs([k - k(w(1)), k - k(w(2))]) + epsilon).^gam;
            Fk(w,:) = ((abs([k - k(w(1)), k - k(w(2))]) + epsilon).^gam)';
        case 'exponential'
            Fk(:,w) = exp((abs([k - k(w(1)), k - k(w(2))]) + epsilon)*gam);
            Fk(w,:) = exp((abs([k - k(w(1)), k - k(w(2))]) + epsilon)*gam)';
        end

    case 'deg-min'

        switch mv2
        case 'powerlaw'
            Fk(:,w) = [min(k,k(w(1))) + epsilon, min(k,k(w(2))) + epsilon].^gam;
            Fk(w,:) = ([min(k,k(w(1))) + epsilon, min(k,k(w(2))) + epsilon].^gam)';
        case 'exponential'
            Fk(:,w) = exp([min(k,k(w(1))) + epsilon, min(k,k(w(2))) + epsilon]*gam);
            Fk(w,:) = exp([min(k,k(w(1))) + epsilon, min(k,k(w(2))) + epsilon]*gam)';
        end

    case 'deg-max'
        
        switch mv2
        case 'powerlaw'
            Fk(:,w) = [max(k,k(w(1))) + epsilon, max(k,k(w(2))) + epsilon].^gam;
            Fk(w,:) = ([max(k,k(w(1))) + epsilon, max(k,k(w(2))) + epsilon].^gam)';
        case 'exponential'
            Fk(:,w) = exp([max(k,k(w(1))) + epsilon, max(k,k(w(2))) + epsilon]*gam);
            Fk(w,:) = exp([max(k,k(w(1))) + epsilon, max(k,k(w(2))) + epsilon]*gam)';
        end

    case 'deg-prod'
        
        switch mv2
        case 'powerlaw'
            Fk(:,w) = ([k*k(w(1)) + epsilon, k*k(w(2)) + epsilon].^gam);
            Fk(w,:) = (([k*k(w(1)) + epsilon, k*k(w(2)) + epsilon].^gam)');
        case 'exponential'
            Fk(:,w) = exp([k*k(w(1)) + epsilon, k*k(w(2)) + epsilon]*gam);
            Fk(w,:) = exp([k*k(w(1)) + epsilon, k*k(w(2)) + epsilon]*gam)';
        end     

    end

    P = Fd.*Fk(indx);
    b(i) = r;
    P(b(1:i)) = 0;
end
b = indx(b);


function b = fcn_nghbrs(A,K,Fd,m,gam,modelvar,epsilon)
K = K + epsilon;
n = length(Fd);
mseed = nnz(A)/2;
A = A > 0;

mv2 = modelvar{2};

switch mv2
    case 'powerlaw'
%         gam = abs(gam);
        Fk = K.^gam;
    case 'exponential'
        Fk = exp(gam*K);
end
Ff = Fd.*Fk.*~A;
[u,v] = find(triu(ones(n),1));
indx = (v - 1)*n + u;
P = Ff(indx);
b = zeros(m,1);
b(1:mseed) = find(A(indx));
for i = (mseed + 1):m
    C = [0; cumsum(P)];
    r = sum(rand*C(end) >= C);
    %save('Test')
    b(i) = r;
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
    switch mv2
        case 'powerlaw'
            Ff(uu,y) = Fd(uu,y).*(K(uu,y).^gam);
            Ff(y,uu) = Ff(uu,y)';
            Ff(vv,x) = Fd(vv,x).*(K(vv,x).^gam);
            Ff(x,vv) = Ff(vv,x)';
        case 'exponential'
            Ff(uu,y) = Fd(uu,y).*exp(K(uu,y)*gam);
            Ff(y,uu) = Ff(uu,y)';
            Ff(vv,x) = Fd(vv,x).*exp(K(vv,x)*gam);
            Ff(x,vv) = Ff(vv,x)';
    end
    Ff(A) = 0;
    P = Ff(indx);
end
b = indx(b);

function b = fcn_matching(A,K,Fd,m,gam,modelvar,epsilon)
K = K + epsilon;
n = length(Fd);
mseed = nnz(A)/2;

mv2 = modelvar{2};

switch mv2
    case 'powerlaw'
        Fk = K.^gam;
    case 'exponential'
        Fk = exp(gam*K);
end
Ff = Fd.*Fk.*~A;
[u,v] = find(triu(ones(n),1));
indx = (v - 1)*n + u;
P = Ff(indx);

b = zeros(m,1);
b(1:mseed) = find(A(indx));

for ii = (mseed + 1):m
    
    C = [0; cumsum(P)];
    r = sum(rand*C(end) >= C);
    b(ii) = r;
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
    switch mv2
        case 'powerlaw'
            Fk = K.^gam;
        case 'exponential'
            Fk = exp(gam*K);
    end
    Ff = Fd.*Fk.*~A;
    P = Ff(indx);
end
b = indx(b);

function b = fcn_sptl(A,Fd,m)
n = length(Fd);
mseed = nnz(A)/2;
[u,v] = find(triu(ones(n),1));
indx = (v - 1)*n + u;
P = Fd(indx).*~A(indx);
b = zeros(m,1);
b(1:mseed) = find(A(indx));
for i = (mseed + 1):m
    C = [0; cumsum(P)];
    r = sum(rand*C(end) >= C);
    b(i) = r;
    P = Fd(indx);
    P(b(1:i)) = 0;
end
b = indx(b);

% What's this? A hidden model not mentioned in the paper!??!? Well yes.
% This does work. In short it calculates the communicability between all 
% pairs of nodes and uses that as the topology value. However, it is very
% VERY slow because every iteration you need to calculate the matrix
% exponential. It didn't make it into the main paper because running it was
% just taking too long and I was on a time limit!
function b = fcn_com(A,K,Fd,m,gam,modelvar,epsilon)
K = K + epsilon;
n = length(Fd);
mseed = nnz(A)/2;
mv2 = modelvar{2};

switch mv2
    case 'powerlaw'
        Fk = K.^gam;
    case 'exponential'
        Fk = exp(K.*gam);
end
Ff = Fd.*Fk.*~A;
[u,v] = find(triu(ones(n),1));
indx = (v - 1)*n + u;
P = Ff(indx);

b = zeros(m,1);
b(1:mseed) = find(A(indx));
for i = (mseed + 1):m
    C = [0; cumsum(P)];
    r = sum(rand*C(end) >= C);
    b(i) = r;
    uu = u(r);
    vv = v(r);
    A(uu,vv) = 1;
    A(vv,uu) = 1;
    switch mv2
        case 'powerlaw'
            Fk = (gexpm(A)+epsilon).^gam;
        case 'exponential' 
            Fk = exp((gexpm(A)+epsilon).*gam);
    end
    Ff = Fd.*Fk.*~A;
    P = Ff(indx);
end
b = indx(b);