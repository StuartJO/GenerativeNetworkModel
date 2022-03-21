function [h,p] = ComputeSigDiff(Vals,sig,BonCorr,CorrType)

% This will compute all the possible correlations between the data in
% "Vals", applying a Bonferroni correciton if desired
%
% Inputs:
%                           Vals = an 1*N cell array where each cell 
%                                  contains data 
%
%                           sig = significance level (scalar)
%
%                           BonCorr = set to 1 to provide a Bonferroni
%                           correlation, 0 to otherwise. To set a custom
%                           number of corrections set BonCorr to an integer
%                           > 1
%
%                           CorrType = 1 for Wilcoxon signed rank test or 2
%                                       for a t-test
%
% Outputs:
%                           h = logical matrix where h(i,j) 1 indicates 
%                           significance (0 otherwise) betweens Vals{i} and 
%                           Vals{j}      
%
%                           p = matrix where p(i,j) indicates the p-value
%                           between betweens Vals{i} and Vals{j}

if nargin < 4
    CorrType = 1;
end

NVals = length(Vals);

if BonCorr == 1
   Ncorr = M*2; 
   Ncorr = Ncorr - Ignr;
elseif BonCorr == 0
    Ncorr = 1;   
else
    Ncorr = BonCorr;
end

p = zeros(NVals);
h = logical(p);

for i = 1:NVals-1
    for j = i+1:NVals
    if isempty(Vals{i}) || isempty(Vals{j})
       p(i,j) = NaN;
       h(i,j) = NaN;
    else
        if CorrType == 1
          [p(i,j),h(i,j)] = signrank(Vals{i},Vals{j},'alpha',sig/Ncorr);
        elseif CorrType == 2
          [h(i,j),p(i,j)] = ttest(Vals{i},Vals{j},'Alpha',sig/Ncorr);  
        end
    end
    end
end
p = p + p';
h = h + h';
     %p = p.*Ncorr;

end