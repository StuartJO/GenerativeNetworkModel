function new_vertices = PointsUnderSurface(surface,d,correct_outside)

% This function finds points located just under the surface at each vertex
% in a triangulated surface. To do this, the vertex normals are calculated
% with a magnitude d. This value is subtracted from each vertex to
% therefore find a point that is d units under the surface, perpendicular
% to the tangent at that vertex

% The "surface" input should be a structure with the fields "vertices" and
% "faces"

% Because the vertex normals can sometimes be wrong (e.g. if a vertex is
% located at a sharp point and there are different numbers of faces on the
% different "sides" of the point, the vertex normal won't be exactly
% perpendicular), an additional check is done. If the new point lies
% outside the surface, we randomly sample points that are d distance away
% from the specific vertex. When one is found that is used as the new
% point.

if nargin < 3
    correct_outside = 0;
end

vertices = surface.vertices;
faces = surface.faces;

nverts = size(vertices,1);

vert1 = vertices(faces(:,1),:);
vert2 = vertices(faces(:,2),:);
vert3 = vertices(faces(:,3),:);
faceNormal = cross(vert2-vert1, vert3-vert1,2);

normals = zeros(size(vertices));

%area = .5*(sqrt(sum(faceNormal.^2, 2)));

for i = 1:nverts
normals(i,:) = sum(faceNormal(logical(sum(double(faces == i),2)),:));
end

mag = sqrt(sum(normals.^2,2));

normals = normals.*(d./[mag mag mag]);

new_vertices = vertices - normals;

if correct_outside
    PT = in_polyhedron(faces, vertices,new_vertices);
    outside = find(PT==0);
    if ~isempty(outside)
        for i = 1:length(outside)
           vert_ind = outside(i);
           vert_x = vertices(vert_ind,1);
           vert_y = vertices(vert_ind,2); 
           vert_z = vertices(vert_ind,3); 
           IN = 0;
           while IN == 0
               r = randn(3,1);
               r = d*bsxfun(@rdivide,r,sqrt(sum(r.^2,1)));
               x = vert_x+r(1,:);
               y = vert_y+r(2,:);
               z = vert_z+r(3,:);
               IN = in_polyhedron(faces, vertices,[x y z]);               
           end
               new_vertices(vert_ind,:) = [x y z];
        end
    end
end