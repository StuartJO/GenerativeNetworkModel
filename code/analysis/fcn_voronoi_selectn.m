function xnew = fcn_voronoi_selectn(x,fval,ndraw,bounds,pow)
% preferentially sample parameter values using voronoi tessellation
%
% Inputs:
%
%   x = [P x n] set of n parameter values
%
%   fval = [P x 1] set of objective function (we treat small values as 
%       indicative of "better" parameters).
%
%   ndraw = number of new samples you want to generate.
%
%   bounds = [P x 2] range of parameters (the bounds in which any 
%       interesting behavior will occur -- can be really broad).
%       bounds(:,1) is the lower limit and bounds(:,2) is the upper limit
%
%   pow = the "alpha" exponent from the paper. If this value is big, then
%       your new samples will be drawn from only the best regions of
%       parameter space. Smaller values will still preferentially sample
%       from good regions of parameter space, but may include some
%       not-so-great samples, as well.
%
% Outputs:
%
%   xnew = [ndraw x n] set of parameters.
%

% we map parameters to the interval [0,1] because there can
% be huge differences of scale

nTotal = size(x,2);

lims = [0,1];

% map parameters
xorig = x;
x = zeros(size(xorig));
IND = 1;
EXCLUDED = [];
INCLUDED = [];
for i = 1:nTotal
	if bounds(i,1) == bounds(i,2)
	    EXCLUDED = [EXCLUDED i];	
	else
	    x(:,IND) = interp1(bounds(i,:),lims,xorig(:,i));
	    IND = IND + 1;
	    INCLUDED = [INCLUDED i];
	end
end

x(:,sum(x)==0) = [];

n = size(x,2);

% perform voronoi tesselation
[v,c] = voronoin(x);

% matlab can generate empty cells, so get rid of those
rmv = cellfun(@isempty,c);
c(rmv) = [];
fval(rmv) = [];

% generate cumulative distribution for each cell based on energy raised to
% negative pow
cdist = [0; cumsum(fval.^-pow)];

% determine which cells you'll draw from -- these will tend to be low
% fval cells if pow is reasonably large.
indpoly = zeros(ndraw,1);
for idraw = 1:ndraw
    indpoly(idraw) = sum(rand*cdist(end) >= cdist);
end
npts = length(fval);
h = hist(indpoly,1:npts);

% loop over the cells
ind = find(h);
count = h(ind);
npoly = length(ind);
xnew = zeros(ndraw,n);
bigcount = 0;

if n == 2

    for i = 1:npoly
    
    % cell index
    indx = c{ind(i)};
    
    % vertices defining cell
    verts = v(indx,:);
    verts(verts(:,1) < lims(1),1) = lims(1);
    verts(verts(:,1) > lims(2),1) = lims(2);
    verts(verts(:,2) < lims(1),2) = lims(1);
    verts(verts(:,2) > lims(2),2) = lims(2);
    
    % get x and y range of cell
    vertsmin = min(verts);
    vertsmax = max(verts);
    
    % generate new samples within this cell
    for j = 1:count(i)
        r = unifrnd(vertsmin(1),vertsmax(1));
        s = unifrnd(vertsmin(2),vertsmax(2));
        check = inpolygon(r,s,verts(:,1),verts(:,2));
        ntries = 0;
        
        % because cells can have weird shapes, check to make sure the
        % points you sample in rectangular space fall within polygon
        while ~check
            r = unifrnd(vertsmin(1),vertsmax(1));
            s = unifrnd(vertsmin(2),vertsmax(2));
            check = inpolygon(r,s,verts(:,1),verts(:,2));
            ntries = ntries + 1;
            if ntries >= (100*count(i))
                break
                % we've never had any issues with this, but if you try too
                % many times, and fail to generate a point within the cell,
                % terminate after some point 
            end
        end
        
        % append new sample to list
        bigcount = bigcount + 1;
        xnew(bigcount,:) = [r,s];
    end
    
    end
    
else

for i = 1:npoly
    
    % cell index
    indx = c{ind(i)};
    
    % vertices defining cell
    verts = v(indx,:);
    for j = 1:n
        verts(verts(:,j) < lims(1),j) = lims(1);
        verts(verts(:,j) > lims(2),j) = lims(2);
    end
    
    % get x and y range of cell
    vertsmin = min(verts);
    vertsmax = max(verts);  
    
    
    % Find the convex hull 
    chull = convhulln(verts);
    
    chull_size = size(chull,1);
    
    cf = zeros(size(chull));
    df = zeros(chull_size,1);
    
    % Find the coefficients of the plane, The planes define the cell of the
    % Voronoi tesselation
    for j = 1:chull_size       
        plane = zeros(n);
        for k = 1:n
            plane(k,:) = verts(chull(j,k),:);
        end
        [cf(j,:),df(j,:)] = plane_coef(plane);
    end
    
    point = zeros(1,n);
    
    % generate new samples within this cell
    for j = 1:count(i)
        
        for k = 1:n
            point(1,k) = unifrnd(vertsmin(k),vertsmax(k));
        end
        
        % This is basically a version of inpolygon but extended to work in
        % n dimensions
        check=~any(cf*point'+df(:,ones(1,size(point,1)))>0,1);
        
        ntries = 0;
        
        % because cells usually have weird shapes, check to make sure the
        % points you sample in rectangular space fall within polygon
        while ~check
            
        for k = 1:n
            point(1,k) = unifrnd(vertsmin(k),vertsmax(k));
        end
        
        check=~any(cf*point'+df(:,ones(1,size(point,1)))>0,1);
        
            ntries = ntries + 1;
            if ntries >= (100*count(i))
                break
                % we've never had any issues with this, but if you try too
                % many times, and fail to generate a point within the cell,
                % terminate after some point
            end
        end
        
        % append new sample to list
        bigcount = bigcount + 1;
        xnew(bigcount,:) = point;
    end
    
end

end

% need to map points back into original parameter space
xneworig = xnew;
xnew = zeros(ndraw,nTotal);

for i = 1:n
    xnew(:,INCLUDED(i)) = interp1(lims,bounds(INCLUDED(i),:),xneworig(:,i));
end

if ~isempty(EXCLUDED)
    for i = 1:length(EXCLUDED)
        xnew(:,EXCLUDED(i)) = bounds(EXCLUDED(i),1);
    end
end


function [cf,df] = plane_coef(plane)
% Calculates the coefficients and constant terms of the equation of a
% plane in implicit form
%
% This function is taken from https://www.mathworks.com/matlabcentral/fileexchange/48509-computational-geometry-toolbox
%

[~,m]=size(plane);
pdiff=diff(plane);
cf=zeros(m,1);
sign=1;
r=1:m;
for i=r
    cf(i)=sign*det(pdiff(:,r~=i));
    sign=-sign;
end
cf=cf/norm(cf);
df=-plane(1,:)*cf;