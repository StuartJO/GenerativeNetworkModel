function [DATA,V,ind] = VisualiseFitsIndv(data1,labels,label_grouping,NBest,ylabelname,plotorder,order,Colors2use,OutlineType)

if nargin < 4
   NBest = []; 
end

if nargin < 5
    ylabelname = 'Model fit';
end

Nmodels = length(data1);

if nargin < 6
    plotorder = 'min';
end

if nargin < 7
    order = 'min';
end

if nargin < 3
    label_grouping = [1 2 2 3 3 3 3 3 4 4 4 4 4 5];
end

for i = 1:length(data1)
if size(data1{i},1) > size(data1{i},2)
    data1{i} = data1{i}';
end
end

if strcmp(ylabelname,'/eta')
    flipSign = -1;
else
    flipSign = 1;
end

if nargin < 8

Colors2use(1,:) = [152,78,163]/255;
Colors2use(2,:) = [55,126,184]/255;
Colors2use(3,:) = [228,26,28]/255;
Colors2use(4,:) = [77,175,74]/255;
Colors2use(5,:) = [255,127,0]/255;
Colors2use(6,:) = [.75,.75,.75];

end
if nargin < 9
    OutlineType = zeros(1,Nmodels);
end

ModelColors = Colors2use(label_grouping,:);

g = cell(1,1);
DATA = cell(1,Nmodels);

switch order
    case 'min'
    NbestOrder ='ascend';
    case 'max'
    NbestOrder ='descend';
    otherwise
     NbestOrder ='ascend';   
end


for i = 1:length(labels)
    
if isempty(NBest)
    
    evals = data1{i};
    MinE(i) = mean(evals);
    DATA{i} = evals;
    g{i} = repmat({[labels{i},' type 1']},length(evals),1);

else
    
    evals = sort(data1{i},NbestOrder);
    
    DATA{i} = evals(1:NBest);
    
    MinE(i) = mean(DATA{i});
    
    g{i} = repmat({join([labels{i},' type 1'])},NBest,1);  

end

end

if isa(plotorder,'double')
    ind = plotorder;
else

    switch plotorder
        case 'min'
        [~,ind] = sort(MinE,'ascend');
        case 'max'
        [~,ind] = sort(MinE,'descend');    
        otherwise
        [~,ind] = sort(MinE,'ascend');    
    end
end

for i = 1:length(ind)
     
   idx(i) = ind(i);
   idx_colors(i) = ind(i);
    
end
    
    DATA_ordered = DATA(idx);
    for i = 1:length(DATA)
    V(1,i) = Violin(DATA_ordered{i}.*flipSign, i,'ShowMean',true);
    V(1,i).ViolinColor = ModelColors(idx_colors(i),:);
    V(1,i).MeanPlot.Color = [0 0 0];
    %if mod(i,2) == 1
    
    %
    OutlineTypeOrdered = OutlineType(idx);
    if OutlineTypeOrdered(i) == 0
    V(1,i).EdgeColor = findAlphaColor(ModelColors(idx_colors(i),:),1);
    elseif OutlineTypeOrdered(i) == 1
        V(1,i).EdgeColor = [0 0 0];
        V(1,i).ViolinPlot.LineWidth = 2;
        V(1,i).ViolinPlot.LineStyle = '-';
    elseif OutlineTypeOrdered(i) == 2
        V(1,i).EdgeColor = [0 0 0];
        V(1,i).ViolinPlot.LineWidth = 2;
        V(1,i).ViolinPlot.LineStyle = '--';
        
    elseif OutlineTypeOrdered(i) == 3
        V(1,i).EdgeColor = [0 0 0];
        V(1,i).ViolinPlot.LineWidth = 2;
        V(1,i).ViolinPlot.LineStyle = ':';
    end
    
        V(1,i).BoxColor = [0 0 0];
    end    
    
    xlim([.5 (Nmodels)+.5])


ylimits = ylim;

xticklabels([]);
xticks(1:Nmodels)
%xticks(1.5:2:((Nmodels*2)+1));

%labels_new_order = labels(ind);

%labels2lines_ordered = cellfun(@(x) strrep(x,' ','\newline'), labels_new_order,'UniformOutput',false);

xticklabels(labels(ind))

%set(gca,'TickLabelInterpreter', 'latex');

%xticklabels(labels2lines_ordered)

xtickangle(45)

ylabel(ylabelname)

xlabel('Model')

ylim(ylimits)

set(gca,'Fontsize',20)
set(gca, 'Layer','top'); set(gca,'TickDir','In')
ylim(ylimits)