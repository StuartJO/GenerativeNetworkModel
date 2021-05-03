function [h,p,val_pairs] = CompareFitsTtest(Vals,sig,BonCorr)

if BonCorr == 1
   N = length(Vals)/2; 
else
    N = 1;   
end

M = length(Vals)/2; 

val_pairs = cell(M,1);

for i = 1:M
    val_pairs{i} = [(i*2)-1 (i*2)];
    [p(i),h(i)] = signrank(Vals{val_pairs{i}(1)},Vals{val_pairs{i}(2)},'alpha',sig/N);
    %[h(i),p(i)] = ttest(Vals{val_pairs{i}(1)},Vals{val_pairs{i}(2)},'Alpha',sig/N);
end

     p = p.*N;


end