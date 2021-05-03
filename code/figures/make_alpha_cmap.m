function cmap = make_alpha_cmap(RGB,N,Scaling)

if nargin < 3

s = linspace(1,0,N)';

else
    s = linspace(1,0,N)';
    s = s.^Scaling;
end

rgb_map = repmat(RGB,N,1);

cmap_adjust = 1-rgb_map;

spacing = [s s s];

cmap = (cmap_adjust.*spacing)+rgb_map;