function [h,p,Ncorr] = CompareFitsSigTest(Vals,sig,BonCorr,CorrType)

if nargin < 4
    CorrType = 1;
end

M = size(Vals,1\2); 

Ignr = 0;

for i = 1:M
    if isempty(Vals{1,i}) || isempty(Vals{2,i})
       Ignr = Ignr + 1; 
    end
end

if BonCorr == 1
   Ncorr = M*2; 
   Ncorr = Ncorr - Ignr;
elseif BonCorr == 0
    Ncorr = 1;   
else
    Ncorr = BonCorr;
end

for i = 1:M
    if isempty(Vals{1,i}) || isempty(Vals{2,i})
       p(i) = 0;
       h(i) = 0;
    else
        if CorrType == 1
          [p(i),h(i)] = signrank(Vals{1,i},Vals{2,i},'alpha',sig/Ncorr);
        elseif CorrType == 2
          [h(i),p(i)] = ttest(Vals{1,i},Vals{2,i},'Alpha',sig/Ncorr);  
        end
    end
end

     %p = p.*Ncorr;

end