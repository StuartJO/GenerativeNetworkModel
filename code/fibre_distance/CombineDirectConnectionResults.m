function [G,Bin,D] = CombineDirectConnectionResults(DATALOC)

%DATALOC = '/projects/kg98/stuarto/Fetal_brain/CRL_Fetal_Brain_Atlas_2017/NewOptimFibreLength';

load([DATALOC,'/Reduced.mat'],'target_points');

N = dlmread([DATALOC,'/ReducedJobs.txt']);

Bin = zeros(size(target_points,1));

D = squareform(pdist(target_points));

for i = 1:N
   DATA = load([DATALOC,'/DirectConnections/Output_',num2str(i),'.mat']); 
   Bin(DATA.source_inds,:) = DATA.DIRECT;      
end

Bin = Bin + Bin';

A = sparse(double(D).*Bin);

Bin = sparse(Bin);

G = graph(full(A));