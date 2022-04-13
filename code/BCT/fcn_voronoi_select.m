function xnew = fcn_voronoi_select(x,e,ndraw,xbounds,ybounds,pow)
% preferentially sample parameter values using voronoi tessellation
%
% Inputs:
%
%   x = [P x 2] set of parameter values (the eta and gamma from the paper).
%   e = [P x 1] set of energies (we treat small values as indicative of
%       "better" samples).
%
%   ndraw = number of new samples you want to generate.
%
%   xbounds = range of eta parameters (the bounds in which any interesting
%       behavior will occur -- can be really broad).
%
%   ybounds = same as above, but for gamma parameters.
%
%   pow = the "alpha" exponent from the paper. If this value is big, then
%       your new samples will be drawn from only the best regions of
%       parameter space. Smaller values will still preferentially sample
%       from good regions of parameter space, but may include some
%       not-so-great samples, as well.
%
% Outputs:
%
%   xnew = [ndraw x 2] set of parameters.
%

% we map parameters to the interval [0,1] because there can
% be huge differences of scale (eta tends to be approx an order of
% magnitude greater than gamma parameters for many models)

xlims = [0,1];
ylims = [0,1];

% map parameters
xorig = x;
x = zeros(size(xorig));
x(:,1) = interp1(xbounds,xlims,xorig(:,1));
x(:,2) = interp1(ybounds,ylims,xorig(:,2));

% perform voronoi tesselation
[v,c] = voronoin(x);

% matlab can generate empty cells, so get rid of those
rmv = cellfun(@isempty,c);
c(rmv) = [];
e(rmv) = [];

% generate cumulative distribution for each cell based on energy raised to
% negative pow
cdist = [0; cumsum(e.^-pow)];

% determine which cells you'll draw from -- these will tend to be low
% energy cells if pow is reasonably large.
indpoly = zeros(ndraw,1);
for idraw = 1:ndraw
    indpoly(idraw) = sum(rand*cdist(end) >= cdist);
end
npts = length(e);
h = hist(indpoly,1:npts);

% loop over the cells
ind = find(h);
count = h(ind);
npoly = length(ind);
xnew = zeros(ndraw,2);
bigcount = 0;
for i = 1:npoly
    
    % cell index
    indx = c{ind(i)};
    
    % vertices defining cell
    verts = v(indx,:);
    verts(verts(:,1) < xlims(1),1) = xlims(1);
    verts(verts(:,1) > xlims(2),1) = xlims(2);
    verts(verts(:,2) < ylims(1),2) = ylims(1);
    verts(verts(:,2) > ylims(2),2) = ylims(2);
    
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
                % terminate after some point (in this 
            end
        end
        
        % append new sample to list
        bigcount = bigcount + 1;
        xnew(bigcount,:) = [r,s];
    end
    
end

% need to map points back into original parameter space
xneworig = xnew;
xnew = zeros(size(xneworig));
xnew(:,1) = interp1(xlims,xbounds,xneworig(:,1));
xnew(:,2) = interp1(ylims,ybounds,xneworig(:,2));