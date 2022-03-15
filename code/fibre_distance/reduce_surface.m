function [surface_out,ind,cdata_out] = reduce_surface(surface,factor,cdata)

if nargin < 3
    cdata = [];
end

% Sometimes the surface variable I use has other fields that don't play
% nice with other MATLAB functions, so we just extract the ones that matter

p.vertices = surface.vertices;
p.faces = surface.faces;

rp = reducepatch(p,factor);

%Create an index of the preserved vertex.
[ind,loc] = ismember(p.vertices,rp.vertices,'rows'); 

%if you want to preserve the index order:
locb = loc(ind);
subind = find(ind);
[~,revsor] = sort(locb);
ind = subind(revsor);

if ~isempty(cdata)
    cdata_out = cdata(ind);
else
    cdata_out = [];
end

surface_out = rp;

end