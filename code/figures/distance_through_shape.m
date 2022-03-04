function [distance, direct_connection, eucl_dist] = distance_through_shape(x,y,quiet,x_shape,y_shape,AllInside)

% This function calcuates the shortest distances between verticies/nodes  
% along the perimeter of a 2D shape. An adjacency matrix is created where
% each element is the Euclidean distance between verticies. Then a line is
% drawn between every pair of points. If that line intersects the perimeter 
% of the shape, those two points are labelled as not having a direct
% connection and in the adjacency matrix the respective element is set to
% 0. When all points have been evaluated, the shortest path between every
% pair of points is then calculated based on the values in the adjacency
% matrix.

if nargin < 3
    quiet = 0;
end

if length(x) ~= length(y)
   error('x and y must be of the same length') 
end

if size(x,1) < size(x,2)
    x = x';  
end

if size(y,1) < size(y,2)
    y = y';  
end

if nargin < 4

    if size(x_shape,1) < size(x_shape,2)
    x_shape = x_shape';  
    end

    if size(y,1) < size(y,2)
        y_shape = y_shape';  
    end
    
    if x_shape(1) ~= x_shape(end)
        x_shape(end+1) = x_shape(1);
    end

    if y_shape(1) ~= y_shape(end)
        y_shape(end+1) = y_shape(1);
    end

    coords = [x_shape(1:end-1) y_shape(1:end-1)];
    
    
else
    
% The coordinates need to form a closed loop so the first pair of
% coordinates need to be repeated
x_shape = x;
y_shape = y;

if x_shape(1) ~= x_shape(end)
    x_shape(end+1) = x_shape(1);
end

if y_shape(1) ~= y_shape(end)
    y_shape(end+1) = y_shape(1);
end

coords = [x_shape(1:end-1) y_shape(1:end-1)];

end

if nargin < 6
    AllInside = 0;
end

% Euclidean distances are first defined

eucl_dist = pdist2(coords,coords);
distance = eucl_dist;

reverseStr = '';

ITER = 1;

% Calculate the midpoints of all coordinates

%[X,Y] = meshgrid(x,y);
[X,Y] = meshgrid(coords(:,1),coords(:,2));

Xt = X';
Yt = Y';

[x_mid,y_mid] = midpoint_coords(X,Yt,Xt,Y);

% Take the midpoint values for the upper triangle and check if these points
% are within the shape

[x_midup,upind] = triu2vec(x_mid,1);
y_midup = triu2vec(y_mid,1);

IN = inpolygon(x_midup,y_midup,x_shape,y_shape);

INDS = find(IN == 1);

% For pairs of coordinates whose midpoints lie outside the shape, set the
% distance to 0

OUT = upind(IN == 0);
distance(OUT) = 0;
unique_vertex_pairs = length(INDS);

% Loop over the remaining coordinates whose midpoints are within the shape 
% and check if these coordinates can be connected entirely within the shape

for i = 1:unique_vertex_pairs
    n = polyxpoly([X(upind(INDS(i))) Xt(upind(INDS(i)))],[Yt(upind(INDS(i))) Y(upind(INDS(i)))], x_shape, y_shape,'unique');
    
    if AllInside
    if length(n) ~= 2
           distance(upind(INDS(i))) = 0;
    end
    else
    if ~isempty(n)
           distance(upind(INDS(i))) = 0;
    end        
    end
    if ~quiet
        msg = sprintf('%d/%d connections evaluated\n', ITER, unique_vertex_pairs);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        ITER = ITER + 1;
    end
end

distance = triu(distance,1) + triu(distance,1)';

direct_connection = double(distance > 0);

distance = graphallshortestpaths(sparse(distance));

end

function [x, y] = midpoint_coords(x1,y1,x2,y2)
    x = (x1+x2)./2;
    y = (y1+y2)./2;
end