function b = fcn_nghbrs(A,K,D,m,eta,gam,modelvar)
K = K + 1e-5;
n = length(D);
mseed = nnz(A)/2;
A = A > 0;
mv1 = modelvar{1};
mv2 = modelvar{2};
switch mv1
    case 'powerlaw'
        Fd = D.^eta;
    case 'exponential'
        Fd = exp(eta*D);
end
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
for i = (mseed + 1):m
    C = [0; cumsum(P)];
    r = sum(rand*C(end) >= C);
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
    %%%%%%%%%%%%%%%%%%%%%
    %     P(isinf(P)) = 0;
    %%%%%%%%%%%%%%%%%%%%%
end
b = find(triu(A,1));