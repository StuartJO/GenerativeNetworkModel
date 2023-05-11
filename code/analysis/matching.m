function m = matching(A)

%A = A>0;

n = length(A);

nei = (A*A).*~eye(n);

deg = sum(A);

% tic
degsum = (deg+deg').*~eye(n);
% toc
% 
% tic
% degmat = repmat(deg,n,1);
% degsum = degmat+degmat';
% toc

m = (nei*2)./( (degsum<=2 & nei~=1) + (degsum-(A.*2))  );