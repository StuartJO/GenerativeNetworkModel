function [maxKS,KS,maxKS_p,KS_p,A_vals,B_vals] = Betzel_energy(A,D,B,customOrMatlab)
% OUTPUTS:
% E: energy

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
