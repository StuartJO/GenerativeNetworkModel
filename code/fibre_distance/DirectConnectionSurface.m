function L = DirectConnectionSurface(surface,targets,sources,type)
% "surface" is a structure with the fields "vertices" and "faces",
% corresponding to an Nv*3 array of Nv vertex coordinates and an Nf*3 array 
% of Nf faces. "targets" is an Nt*3 array of Nt target coordinates and
% "sources" is an Ns*3 array of coordinates to measure the fibre distances
% to the target coordinates. "type" indicates what method to use to
% calculate fibre distances. type == 1 is a slower method which assumes the
% target/source coordinates overlap with vertex coordinates. It does dome
% extra steps to ensure the ray intersection doesn't run into precision
% errors. type == 2 assumes the target/source coordinates do not overlap
% with the vertex coordinates. This is faster.

vertices = surface.vertices;
faces = surface.faces;

if nargin < 4
    type = 1;
end


RMvertexfaces = 1;

Ntargets = size(targets,1);
Nsources = size(sources,1);

L = zeros(Nsources,Ntargets);

%normals = vertexNormal(vertices, faces);

vert1 = vertices(faces(:,1),:);
vert2 = vertices(faces(:,2),:);
vert3 = vertices(faces(:,3),:);

Nfaces = size(faces,1);

if type == 1
% type = 1 assumes the points are the vertices themselves
for S = 1:Nsources
    tic
    orig = [sources(S,1) sources(S,2) sources(S,3)];
        
    for T = 1:Ntargets
    %tic

    targ = [targets(T,1) targets(T,2) targets(T,3)];

    dir = targ - orig;

    [intersect] = TriangleRayIntersection (orig, dir, vert1, vert2, vert3,'lineType','ray','planeType','two sided');

    [intersect1] = TriangleRayIntersection (targ, dir, vert1, vert2, vert3,'lineType','ray','planeType','two sided');

    intersects_total = double(intersect-intersect1);
    
    %This check if an intersection is detected along any of the faces to
    %which the source and target vertices belong. If one is detected this
    %is likely a precision error and so will ignore that intersection
    if RMvertexfaces
        
        Sfaces = logical(sum(double(faces == S),2));
        Tfaces = logical(sum(double(faces == T),2));

        intersects_total(Sfaces) = 0;
        intersects_total(Tfaces) = 0;
    end
    % Check how many intersections were detected
    indirectcon = sum(intersects_total);
    if indirectcon == 0
        if RMvertexfaces
            if max(Sfaces+Tfaces) > 1
                directcon = 1;            
            else
                midpoint = (orig + targ)/2;
                directcon = in_polyhedron(faces, vertices,midpoint);
            end    
        else
            midpoint = (orig + targ)/2;
            directcon = in_polyhedron(faces, vertices,midpoint);        
        end
        if directcon == 1
            L(S,T) = 1;
        end
    end
    end
    toc
end

else

for S = 1:Nsources
    tic
       
    orig = repmat(sources(S,:),Nfaces,1);
    for T = 1:Ntargets
        targ = [targets(T,1) targets(T,2) targets(T,3)];
        dir = targ - orig;
    
    intersect = TriangleRayIntersection (orig, dir, vert1, vert2, vert3,'lineType','segment','planeType','two sided');
        if sum(intersect) == 0
            L(S,T) = 1;
        end
    end
    toc
end    
    
end