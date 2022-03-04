function [avgCDFs,CDFvals] = mean_ecdf(DATA)

if iscell(DATA)
    for i = 1:length(DATA)
        
        minVals(i) = min(DATA{i}(:));
        maxVals(i) = max(DATA{i}(:));
        
    end
    
    minVal = min(minVals);
    maxVal = max(maxVals);
    
    n = length(DATA);

    CDFvals = linspace(minVal,maxVal,1000);

avgCDFs = zeros(size(CDFvals));
ttlCDFs = zeros(n,length(avgCDFs));

for i = 1:n

for iCDF = 1:numel(avgCDFs)  % Probably can be vectorized for speed, but this shows the idea
   ttlCDFs(i,iCDF) = mean((DATA{i}<=CDFvals(iCDF)));
end

end
    
else
minVal = min(DATA(:));
maxVal = max(DATA(:));

n = size(DATA,1);
CDFvals = linspace(minVal,maxVal,1000);

avgCDFs = zeros(size(CDFvals));
ttlCDFs = zeros(n,length(avgCDFs));

for i = 1:n

for iCDF = 1:numel(avgCDFs)  % Probably can be vectorized for speed, but this shows the idea
   ttlCDFs(i,iCDF) = mean((DATA(i,:)<=CDFvals(iCDF)));
end

end

end

avgCDFs = mean(ttlCDFs);  % Same as the number of passes through the loop
%plot(CDFvals,avgCDFs)