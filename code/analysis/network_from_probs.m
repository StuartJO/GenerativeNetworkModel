function adj = network_from_probs(prob,desired_density,existing)
     
%   if eta < 0
%     warning('Should specify positive eta?!')
%   end
%   eta = -eta;

  NODES = length(prob);
  %initialize matrices of adjacency, common degrees, and probabilitiy of links

  new_connections = zeros(NODES,NODES);
  
  prob(1:NODES+1:NODES*NODES)=0;   % remove infinities from diagonal (and set it to zero)
  if exist('existing','var')
      prob = prob.*(1-existing);
      adj = existing;
  else
      adj = zeros(NODES,NODES);   % connectivity matrix (no distances!) 
  end
  
if desired_density < 1

  desiredEdges = ceil(desired_density*NODES*(NODES-1)/2);
else

  desiredEdges = desired_density;

end

  ExistingEdges = nnz(triu(adj));
  EdgesToMake = desiredEdges - ExistingEdges;
  
  % Get a vector of each edges position in the matrix and it's associated
  % probability (upper triangle only)
  
  [probvec,ind] = triu2vec(prob,1);
  
  % Randomly sample edges without replacement. Edges are weighted by
  % probvec./sum(probvec). Edges that already exist have been set to have a
  % probability of zero thus these edges can never be chosen

  %NewEdgesIndTriu=randsampleWRW(ind,EdgesToMake,probvec);
  NewEdgesIndTriu = datasample(ind,EdgesToMake,'Replace',false,'Weights',probvec);

% A more complex but slightly slower way to make the matrix. As only the 
% upper triangle is used in the above calculation, only the index of new 
% connections in the upper triangle is returned. The below code finds the
% corresponding index of lower triangles connections

%   NODESADDONE = NODES+1;
%   distfromdiag = 1+ abs(NewEdgesIndTriu - NODESADDONE.*ceil(NewEdgesIndTriu./NODESADDONE));
%   NewEdgesIndTril = NewEdgesIndTriu - distfromdiag.*(NODES-1);
%   adj([NewEdgesIndTriu NewEdgesIndTril]) = 1;

  new_connections(NewEdgesIndTriu) = 1;
  new_connections = new_connections + new_connections';
  
  adj = adj + new_connections;

% Old code, slower than new implementation  
%   E = nnz(triu(adj));
%   %start adding edges one by one, upto a density of links = "Cost"
%   for Enum=E+1:ceil(desired_density*NODES*(NODES-1)/2)
%       
%        %pick an edge weighted by the probabilities of all edges
%        [Erow, Ecol, plist]=find(prob);      %find all the non-zero probabilities
%        ind=1:length(plist);
%        Edges=randsampleWRW(ind,1,plist);    %sample one edge at random, weighted by the probabilities
%        %and set that edge to '1' in the adjacency matrix
%        n=Erow(Edges);
%        m=Ecol(Edges);
%        adj(n,m)=1;
%        adj(m,n)=1; 
%                      
%             %disallow self-links
%             for i=1:NODES
%             prob(i,i)=0;
%             end
%             %set to zero the probability of picking the same edges again later
%             prob(adj==1) = 0;
%               
%   end
end

function [vec,ind] = triu2vec(mat,k)

% This function simply gets the upper triangle of a matrix and converts it
% to a vector
onesmat = ones(size(mat));
UT = triu(onesmat,k);
vec = mat(UT == 1);
ind = find(UT == 1);
end