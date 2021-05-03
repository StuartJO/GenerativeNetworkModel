function [E,KSstat,Ep,KSpVal] = Betzel_energy(A,D,Asim,customOrMatlab)
% OUTPUTS:
% E: energy

if size(A,1) ~= size(Asim,1) || size(A,2) ~= size(Asim,2)
   Asim_mat = zeros(size(A));
   Asim_mat(Asim) = 1;
   Asim_mat = Asim_mat + Asim_mat';
   Asim = Asim_mat;
end

if nargin < 4
    customOrMatlab = 'matlab';
end
%-------------------------------------------------------------------------------

x = cell(4,1);
x{1} = sum(A,2);
x{2} = clustering_coef_bu(A);
x{3} = betweenness_bin(A)';
x{4} = D(triu(A,1) > 0);

y = cell(4,1);
y{1} = sum(Asim,2);
y{2} = clustering_coef_bu(Asim);
y{3} = betweenness_bin(Asim)';
y{4} = D(triu(Asim,1) > 0);

% Compute ks-statistic for each network property -> K
numStats = 4;
KSstat = zeros(1,numStats);
KSpVal = zeros(1,numStats);
for j = 1:numStats
    switch customOrMatlab
    case 'matlab'
        [~,KSpVal(j),KSstat(j)] = kstest2(x{j},y{j});
    case 'custom'
        KSstat(j) = fcn_ks(x{j},y{j});
        KSpVal(j) = 0;
    end
end

[E,idx] = max(KSstat);
Ep = KSpVal(idx);

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
