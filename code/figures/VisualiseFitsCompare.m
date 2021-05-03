function [out_data,V,ind] = VisualiseFitsCompare(data1,data2,labels,label_grouping,ylabelname,order,Colors2use)

%figure('Position',[229 416 1370 688]);

if nargin < 5
    ylabelname = 'Model fit';
end

Nmodels = length(data1);

if nargin < 6
    order = 'min';
end

if nargin < 4
    label_grouping = [1 2 2 3 3 3 3 3 4 4 4 4 4 5];
end

for i = 1:length(data1)
if size(data1{i},1) > size(data1{i},2)
    data1{i} = data1{i}';
end
if size(data2{i},1) > size(data2{i},2)
    data2{i} = data2{i}';
end
end

if strcmp(ylabelname,'\eta')
    flipSign = -1;
else
    flipSign = 1;
end

if nargin < 7

Colors2use(1,:) = [152,78,163]/255;
Colors2use(2,:) = [55,126,184]/255;
Colors2use(3,:) = [228,26,28]/255;
Colors2use(4,:) = [77,175,74]/255;
Colors2use(5,:) = [255,127,0]/255;
%Colors2use(6,:) = [.5,.5,.5];
Colors2use(6,:) = [177,89,40]/255;
Colors2use(7,:) = [255,172,228]/255;
end

ModelColors = Colors2use(label_grouping,:);

DATA = cell(1,Nmodels);
for i = 1:length(labels)
    

    
    evals = data1{i};
    MinE(i) = mean(evals);
    DATA{(i*2)-1} = evals;
    
    evals = data2{i};
    DATA{(i*2)} = evals;  
    


end

if isa(order,'double')
    ind = order;
else

switch order
    case 'min'
    [~,ind] = sort(MinE,'ascend');
    case 'max'
    [~,ind] = sort(MinE,'descend');    
    otherwise
    [~,ind] = sort(MinE,'ascend');    
end

end

for i = 1:length(ind)
     
   idx(((i*2)-1):i*2) = ((ind(i)*2)-1):ind(i)*2;
   idx_colors(((i*2)-1):i*2) = ind(i);
    
end

    out_data = DATA(idx);

    VInd = 1;
    for i = 1:length(DATA)
        
        if ~isempty(out_data{i})
        
            if mod(i,2) == 1
                if isempty(out_data{i+1})
                    centre_violin = .5;
                else
                   centre_violin = 0; 
                end
            else
                 centre_violin = 0;
            end
            
    V(1,VInd) = Violin(out_data{i}.*flipSign, i+centre_violin,'ShowMean',true,'ShowData',false,'ViolinAlpha',1);
    V(1,VInd).ViolinColor = ModelColors(idx_colors(i),:);
    V(1,VInd).MeanPlot.Color = [0 0 0];
    if mod(i,2) == 1
    %V(1,i).EdgeColor = findAlphaColor(ModelColors(idx_colors(i),:),1);
        V(1,VInd).EdgeColor = [0 0 0];
        V(1,VInd).ViolinPlot.LineWidth = 2;
        V(1,VInd).ViolinPlot.LineStyle = '-';    
    else
        V(1,VInd).EdgeColor = [0 0 0];
        V(1,VInd).ViolinPlot.LineWidth = 2;
        V(1,VInd).ViolinPlot.LineStyle = ':';
    end
    V(1,VInd).BoxColor = [0 0 0];
    VInd = VInd+1;
        else
    
        end
    
    end    

    xlim([.5 (Nmodels*2)+.5])


ylimits = ylim;

xticklabels([]);
xticks(1.5:2:((Nmodels*2)+1));

%labels_new_order = labels(ind);

%labels2lines_ordered = cellfun(@(x) strrep(x,' ','\newline'), labels_new_order,'UniformOutput',false);

xticklabels(labels(ind))

%set(gca,'TickLabelInterpreter', 'latex');

%xticklabels(labels2lines_ordered)

xtickangle(45)

ylabel(ylabelname)

xlabel('Topological term')

ylim(ylimits)

for i = 1:length(ind)
    plot([-.5 -.5]+((2*ind(i))-1),[-1 1],'Color','k','LineWidth',1.5);
    plot([.5 .5]+(2*ind(i)),[-1 1],'Color','k','LineWidth',1.5);  
end

set(gca,'Fontsize',20)
set(gca, 'Layer','top'); set(gca,'TickDir','In')
ylim(ylimits)
%LegendPos = L.Position;
%L.ItemTokenSize(2) = 10;