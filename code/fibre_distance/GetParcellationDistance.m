function dists = GetParcellationDistance(DATADIR,parcdata)

    A = zeros(max(parcdata));
    load([DATADIR,'/Reduced.mat'],'vertex_ind')
    
    D = zeros(length(ind));
    FILES = dir([DATADIR,'/Output*']);
    
    for j = 1:length(FILES)
       DistMatPart = load([DATADIR,'/Output_',num2str(j),'.mat']);
       D(DistMatPart.NODES,:) = DistMatPart.DIST;  
    end
    
    PARC_INDS = parcdata(vertex_ind);
    % For each pair of regions, loop over all the distances between
    % vertices which make up those respective ROIs and calculate the mean
    for j = 1:max(parcdata)-1
            JINDS = PARC_INDS == j;
        for k = j+1:max(parcdata)
            KINDS = PARC_INDS == k;
            
            distsjk = D(JINDS,KINDS);
            
            A(j,k) = mean(distsjk(:));
                       
        end
        
    end
    
    dists = A + A';
    

end