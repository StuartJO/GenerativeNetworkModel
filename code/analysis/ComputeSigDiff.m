function [h,p] = ComputeSigDiff(Vals,sig,BonCorr,CorrType)

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