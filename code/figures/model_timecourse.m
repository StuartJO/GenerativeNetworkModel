function [FA,D,K,FD,FK] = model_timecourse(b,Dists,modeltype,modelvar,eta,gam,epsilon,seedNet,normtype)

    n = length(Dists);

    if isempty(seedNet)
        
        A = zeros(n);
        
    else
    
        A = seedNet;
    
    end
    
    mseed = nnz(A)/2;
    Edges2Add = b(mseed+1:end);
    m = length(Edges2Add);
    mv1 = modelvar{1};
    mv2 = modelvar{2};
    
    [u,v] = find(triu(ones(n),1));
    indx = (v - 1)*n + u;

    FA = zeros(m,length(indx));
    D = FA;
    K = FA;
    FD = FA;
    FK = FA;
    
switch modeltype
    
    case 'clu-avg'
        
        for i = 1:m
            
        clu = clustering_coef_bu(A);
        FA(i,:) = A(indx);
        Kseed = (bsxfun(@plus,clu(:,ones(1,n)),clu')/2)+epsilon;
        K(i,:) = Kseed(indx);
        fk = norm_term(Kseed,gam,mv2,normtype(2),A);
        FK(i,:) = fk(indx);
        D(i,:) = Dists(indx);
        fd = norm_term(Dists,eta,mv1,normtype(1),A);
        FD(i,:) = fd(indx);
        
        A = zeros(n);
        A(b(1:mseed+i)) = 1;
        A = A + A';
                       
        end

    case 'clu-diff'

        for i = 1:m
            
        clu = clustering_coef_bu(A);
        FA(i,:) = A(indx);
        Kseed = abs(bsxfun(@minus,clu(:,ones(1,n)),clu'))+epsilon;
        K(i,:) = Kseed(indx);
        fk = norm_term(Kseed,gam,mv2,normtype(2),A);
        FK(i,:) = fk(indx);
        D(i,:) = Dists(indx);
        fd = norm_term(Dists,eta,mv1,normtype(1),A);
        FD(i,:) = fd(indx);

        A = zeros(n);
        A(b(1:mseed+i)) = 1;
        A = A + A';
        
        end

    case 'clu-max'

        for i = 1:m
            
        clu = clustering_coef_bu(A);
        FA(i,:) = A(indx);
        Kseed = bsxfun(@max,clu(:,ones(1,n)),clu')+epsilon;
        K(i,:) = Kseed(indx);
        fk = norm_term(Kseed,gam,mv2,normtype(2),A);
        FK(i,:) = fk(indx);
        D(i,:) = Dists(indx);
        fd = norm_term(Dists,eta,mv1,normtype(1),A);
        FD(i,:) = fd(indx);
        
        A = zeros(n);
        A(b(1:mseed+i)) = 1;
        A = A + A';
                       
        end

    case 'clu-min'

        for i = 1:m
            
        clu = clustering_coef_bu(A);
        FA(i,:) = A(indx);
        Kseed = bsxfun(@min,clu(:,ones(1,n)),clu')+epsilon;
        K(i,:) = Kseed(indx);
        fk = norm_term(Kseed,gam,mv2,normtype(2),A);
        FK(i,:) = fk(indx);
        D(i,:) = Dists(indx);
        fd = norm_term(Dists,eta,mv1,normtype(1),A);
        FD(i,:) = fd(indx);
        
        A = zeros(n);
        A(b(1:mseed+i)) = 1;
        A = A + A';
                       
        end

    case 'clu-prod'

        for i = 1:m
            
        clu = clustering_coef_bu(A);
        FA(i,:) = A(indx);
        Kseed = (clu*clu')+epsilon;
        K(i,:) = Kseed(indx);
        fk = norm_term(Kseed,gam,mv2,normtype(2),A);
        FK(i,:) = fk(indx);
        D(i,:) = Dists(indx);
        fd = norm_term(Dists,eta,mv1,normtype(1),A);
        FD(i,:) = fd(indx);
        
        A = zeros(n);
        A(b(1:mseed+i)) = 1;
        A = A + A';
                       
        end

    case 'deg-avg'

        for i = 1:m

        FA(i,:) = A(indx);
        kseed = sum(A,2);
        Kseed = (bsxfun(@plus,kseed(:,ones(1,n)),kseed')/2)+epsilon;
        K(i,:) = Kseed(indx);
        fk = norm_term(Kseed,gam,mv2,normtype(2),A);
        FK(i,:) = fk(indx);
        D(i,:) = Dists(indx);
        fd = norm_term(Dists,eta,mv1,normtype(1),A);
        FD(i,:) = fd(indx);
        
        A = zeros(n);
        A(b(1:mseed+i)) = 1;
        A = A + A';
                       
        end

    case 'deg-diff'

        for i = 1:m
            
        FA(i,:) = A(indx);
        kseed = sum(A,2);
        Kseed = abs(bsxfun(@minus,kseed(:,ones(1,n)),kseed'))+epsilon;
        K(i,:) = Kseed(indx);
        fk = norm_term(Kseed,gam,mv2,normtype(2),A);
        FK(i,:) = fk(indx);
        D(i,:) = Dists(indx);
        fd = norm_term(Dists,eta,mv1,normtype(1),A);
        FD(i,:) = fd(indx);
        
        A = zeros(n);
        A(b(1:mseed+i)) = 1;
        A = A + A';
                       
        end

    case 'deg-max'

        for i = 1:m
            
        FA(i,:) = A(indx);
        kseed = sum(A,2);
        Kseed = (bsxfun(@max,kseed(:,ones(1,n)),kseed'))+epsilon;
        K(i,:) = Kseed(indx);
        fk = norm_term(Kseed,gam,mv2,normtype(2),A);
        FK(i,:) = fk(indx);
        D(i,:) = Dists(indx);
        fd = norm_term(Dists,eta,mv1,normtype(1),A);
        FD(i,:) = fd(indx);
        
        A = zeros(n);
        A(b(1:mseed+i)) = 1;
        A = A + A';
                       
        end

    case 'deg-min'

        for i = 1:m
            
        FA(i,:) = A(indx);
        kseed = sum(A,2);
        Kseed = (bsxfun(@min,kseed(:,ones(1,n)),kseed'))+epsilon;
        K(i,:) = Kseed(indx);
        fk = norm_term(Kseed,gam,mv2,normtype(2),A);
        FK(i,:) = fk(indx);
        D(i,:) = Dists(indx);
        fd = norm_term(Dists,eta,mv1,normtype(1),A);
        FD(i,:) = fd(indx);
        
        A = zeros(n);
        A(b(1:mseed+i)) = 1;
        A = A + A';
                       
        end

    case 'deg-prod'

        for i = 1:m

        FA(i,:) = A(indx);
        kseed = sum(A,2);
        Kseed = ((kseed*kseed').*~eye(n))+epsilon;
        K(i,:) = Kseed(indx);
        fk = norm_term(Kseed,gam,mv2,normtype(2),A);
        FK(i,:) = fk(indx);
        D(i,:) = Dists(indx);
        fd = norm_term(Dists,eta,mv1,normtype(1),A);
        FD(i,:) = fd(indx);
        
        A = zeros(n);
        A(b(1:mseed+i)) = 1;
        A = A + A';
                       
        end

    case 'neighbors'

        for i = 1:m
            
        FA(i,:) = A(indx);
        Kseed = ((A*A).*~eye(n))+epsilon;
        K(i,:) = Kseed(indx);
        fk = norm_term(Kseed,gam,mv2,normtype(2),A);
        FK(i,:) = fk(indx);
        D(i,:) = Dists(indx);
        fd = norm_term(Dists,eta,mv1,normtype(1),A);
        FD(i,:) = fd(indx);
        
        A = zeros(n);
        A(b(1:mseed+i)) = 1;
        A = A + A';
                       
        end

    case 'matching'

        for i = 1:m
            
        FA(i,:) = A(indx);
        Kseed = matching_ind(A);
        Kseed = (Kseed + Kseed')+epsilon;
        K(i,:) = Kseed(indx);
        fk = norm_term(Kseed,gam,mv2,normtype(2),A);
        FK(i,:) = fk(indx);
        D(i,:) = Dists(indx);
        fd = norm_term(Dists,eta,mv1,normtype(1),A);
        FD(i,:) = fd(indx);
        
        A = zeros(n);
        A(b(1:mseed+i)) = 1;
        A = A + A';
                       
        end

    case 'sptl'
        for i = 1:m
            
        FA(i,:) = A(indx);
        D(i,:) = Dists(indx);
        fd = norm_term(Dists,eta,mv1,normtype(1),A);
        FD(i,:) = fd(indx);
        
        A = zeros(n);
        A(b(1:mseed+i)) = 1;
        A = A + A';
                       
        end

    case 'com'

        for i = 1:m
            
        FA(i,:) = A(indx);
        Kseed = gexpm(A)+epsilon;
        K(i,:) = Kseed(indx);
        fk = norm_term(Kseed,gam,mv2,normtype(2),A);
        FK(i,:) = fk(indx);
        D(i,:) = Dists(indx);
        fd = norm_term(Dists,eta,mv1,normtype(1),A);
        FD(i,:) = fd(indx);
        
        A = zeros(n);
        A(b(1:mseed+i)) = 1;
        A = A + A';
                       
        end

end