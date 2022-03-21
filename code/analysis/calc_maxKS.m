function [maxKS,KS,maxKS_p,KS_p,A_vals,B_vals] = calc_maxKS(A,D,B,customOrMatlab)

% This function calculates the maxKS value for network B when comparing it
% to network A
%
% Inputs:
%                       A = an N*N adjacency matrix
%
%                       D = an N*N distance matrix (should correspond to
%                       distances between nodes in A/B)
%           
%                       B = an N*N adjacency matrix to compare to A. B can
%                       also be an list of edge indices
%
%                       customOrMatlab = 'matlab' (default) or 'custom'.
%                       Uses either MATLABs inbult function to calculate KS
%                       values or one provided in this script (both will
%                       give the same result, however the MATLAB one will
%                       provide p-values)
% Outputs:
%                       maxKS = the maximum KS value (i.e., indication fo
%                       model fit)
%
%                       KS = KS values for each measure in a 1*4 array,
%                       where KS(1) is the KS value for the degree
%                       distribution, KS(2) is the clustering distribution,
%                       KS(3) is the betweenness distribution, KS(4) is the
%                       edge length distribution 
%
%                       maxKS_p = p-values for each KS value (only if 
%                       customOrMatlab = 'matlab')
%                       
%                       A_vals = the values for input network A
%
%                       B_vals = the values for input network B

if size(D,1) ~= size(B,1) || size(D,2) ~= size(B,2)
   Asim_mat = zeros(size(D));
   Asim_mat(B) = 1;
   Asim_mat = Asim_mat + Asim_mat';
   B = Asim_mat;
end

if nargin < 4
    customOrMatlab = 'matlab';
end
%-------------------------------------------------------------------------------
if ~iscell(A)
A_vals = cell(4,1);
A_vals{1} = sum(A,2);
A_vals{2} = clustering_coef_bu(A);
A_vals{3} = betweenness_bin(A)';
A_vals{4} = D(triu(A,1) > 0);
else
   A_vals = A; 
end

B_vals = cell(4,1);
B_vals{1} = sum(B,2);
B_vals{2} = clustering_coef_bu(B);
B_vals{3} = betweenness_bin(B)';
B_vals{4} = D(triu(B,1) > 0);

% Compute ks-statistic for each network property -> K
numStats = 4;
KS = zeros(1,numStats);
KS_p = zeros(1,numStats);
for j = 1:numStats
    switch customOrMatlab
    case 'matlab'
        [~,KS_p(j),KS(j)] = kstest2(A_vals{j},B_vals{j});
    case 'custom'
        KS(j) = fcn_ks(A_vals{j},B_vals{j});
        KS_p(j) = 0;
    end
end

[maxKS,idx] = max(KS);
maxKS_p = KS_p(idx);

%-------------------------------------------------------------------------------
function kstat = fcn_ks(x1,x2)
    binEdges    =  [-inf ; sort([x1;x2]) ; inf];

    binCounts1  =  histc(x1, binEdges, 1);
    binCounts2  =  histc(x2, binEdges, 1);

    sumCounts1  =  cumsum(binCounts1)./sum(binCounts1);
    sumCounts2  =  cumsum(binCounts2)./sum(binCounts2);

    sampleCDF1  =  sumCounts1(1:end-1);
    sampleCDF2  =  sumCounts2(1:end-1);

    deltaCDF  =  abs(sampleCDF1 - sampleCDF2);
    kstat = max(deltaCDF);
end
%-------------------------------------------------------------------------------
end
