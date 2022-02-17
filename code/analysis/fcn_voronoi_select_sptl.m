function xnew = fcn_voronoi_select_sptl(x,e,ndraw,xbounds,pow)

xlims = [0,1];
% ylims = [0,1];

xorig = x;
x = interp1(xbounds,xlims,xorig(:,1));
[xsort,indsort] = sort(x,'ascend');
esort = e(indsort);

bounds = [0; (xsort(2:end) + xsort(1:end - 1))/2; 1];

cdist = [0; cumsum(esort.^-pow)];
indpoly = zeros(ndraw,1);
for idraw = 1:ndraw
    indpoly(idraw) = sum(rand*cdist(end) >= cdist);
end
npts = length(e);
h = hist(indpoly,1:npts);

ind = find(h);
count = h(ind);
npoly = length(ind);
xnew = zeros(ndraw,1);
bigcount = 0;
for i = 1:npoly
    
    verts = bounds(ind(i):(ind(i) + 1));
    
    vertsmin = min(verts);
    vertsmax = max(verts);
    
    r = unifrnd(vertsmin(1),vertsmax(1),count(i),1);
    xnew((bigcount + 1):(bigcount + count(i))) = r;
    bigcount = bigcount + count(i);
%     
%     for j = 1:count(i)
%         r = unifrnd(vertsmin(1),vertsmax(1));
%         s = unifrnd(vertsmin(2),vertsmax(2));
%         check = inpolygon(r,s,verts(:,1),verts(:,2));
%         while ~check
%             r = unifrnd(vertsmin(1),vertsmax(1));
%             s = unifrnd(vertsmin(2),vertsmax(2));
%             check = inpolygon(r,s,verts(:,1),verts(:,2));
%         end
%         bigcount = bigcount + 1;
%         xnew(bigcount,:) = [r,s];
%     end
    
end
xneworig = xnew;
% xnew = zeros(size(xneworig));
xnew = interp1(xlims,xbounds,xneworig(:,1));
% xnew(:,2) = interp1(ylims,ybounds,xneworig(:,2));

