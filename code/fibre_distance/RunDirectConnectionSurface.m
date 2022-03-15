function [DIRECT,source_inds] = RunDirectConnectionSurface(surface,target_points,source_inds)

% surface = structure with fields "vertices" and "faces"
% target_points = an N*3 list of coordinates to check for direct
% connections
% source_inds = an array of indices indicating which target_points to check
% against all others. Needs to be a consecutive list of integers

Ntarget_points = size(target_points,1);

IsConsecutive = unique(diff(source_inds));

if length(IsConsecutive) ~= 1 && IsConsecutive(1) == 1
   error('''source_inds'' needs to be consecutive integers!') 
end

% If the number of points requested exceeds the number of vertices, only
% compute for the vertices that exist!
source_inds(source_inds > Ntarget_points-1)=[];

Nsource_points = length(source_inds);

DIRECT=zeros(Nsource_points,Ntarget_points); 
for i = 1:Nsource_points
    source_ind = source_inds(i);
    source_point = target_points(source_ind,:);
    target_inds = source_ind+1:Ntarget_points;
    DIRECT(i,target_inds) = DirectConnectionSurface(surface,target_points(target_inds,:),source_point,2);
end 