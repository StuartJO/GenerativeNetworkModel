function [ordered_mdl_data,V,mdl_order,issig,pval] = PlotMdlResults(MdlData,MdlNames,varargin)

%[ordered_mdl_data,V,mdl_order,h,p] = PlotPairedMdlResults(DataMdlGrp1,DataMdlGrp2,MdlNames,GrpNames,DataLabel,MdlOrder,MdlTypesInd,MdlTypesCmap,MdlTypesNames,SigLvl)

% This function allows you to compare pairs of distributions for different
% models.
%
% INPUTS:
%
% mdl_grp1_data = a cell array (1*N, where N is the number of models).
% Each cell should contain a vector of model outputs (e.g., maxKS, Fcv,
% eta/gam values etc) for a set of grouped models
%
% mdl_grp2_data = same format as mdl_grp1_data. Each cell in it should
% correspond to the same model type as that in mdl_pair1_data
%
% mdl_names = a cell array with the names of each of the N models
%
% mdl_subgrp = a vector of size 1*N indicatign for each model to which sub
% grouping it belongs
%
% ylabelname = a string indicating the name of what the data being plotted
% is
%
% order = indicates the order in which model's are to be plotted. Can
% either be 'min' (model with the smallest mean value is plotted first) or
% 'max' (model with the largest mean value is plotted first), or a 1*N
% vector indicating the order (starting from 1) with which to plot each
% paired set of models
%
% mdl_subgrp_colors = a colormap indicating which colours to associated
% with each sub group (mdl_subgrp_colors(1,:) will be the color associated
% with mdl_subgrp == 1 etc etc)
%
% mdl_subgrp_names = a cell containing the name of each model sub group
%
% mdl_grp_names = a cell containing the name of the two model groupings
%
% SigLvl = indicates if pairs of grouped models should be assessed for
% statistical significance and significant plotted, and how that
% significance is determined. If SigLvl = 1, significance will be a
% Bonferroni corrected Wilcoxon signed rank test, controlling for N-n tests
% (where n is the number of paired models for which only one set of data is
% being plotted). If SigLvl < 1, an uncorrected Wilcoxon signed rank test
% with a p value = SigLvl will be used. If SigLvl > 1, a Bonferroni 
% corrected Wilcoxon signed rank test controlling for SigLvl tests is
% performed. If SigLvl = 0, no significance testing will be performed.

Nmodels = size(MdlData,2);
Ngrps = size(MdlData,1);

defaultCmap = [152,78,163;...
55,126,184;...
228,26,28;...
77,175,74;...
255,127,0;...
66,224,245;...
166,86,40;...
77,190,238]./255;

p = inputParser;

defaultDataLabel = '';
defaultMdlOrder = 'min';
defaultSigLevel = 0;

defaultGrpNames = {'Model group A','Model group B'};

defaultMdlTypesInd = ones(1,Nmodels);
defaultMdlTypesCmap = defaultCmap;
defaultMdlTypesNames = [];

addRequired(p,'MdlData',@(x) iscell(x))
addRequired(p,'MdlNames',@(x) iscell(x))
addParameter(p,'GrpNames',defaultGrpNames,@(x) iscell(x))
addParameter(p,'DataLabel',defaultDataLabel,@(x) ischar(x))
addParameter(p,'MdlTypesInd',defaultMdlTypesInd)
addParameter(p,'MdlTypesCmap',defaultMdlTypesCmap)
addParameter(p,'MdlTypesNames',defaultMdlTypesNames)
addParameter(p,'MdlOrder',defaultMdlOrder)
addParameter(p,'SigLvl',defaultSigLevel)

parse(p,MdlData,MdlNames,varargin{:});

MdlData = p.Results.MdlData;
MdlNames = p.Results.MdlNames;
GrpNames = p.Results.GrpNames;
DataLabel = p.Results.DataLabel;
MdlTypesInd = p.Results.MdlTypesInd;
MdlTypesCmap = p.Results.MdlTypesCmap;
MdlTypesNames = p.Results.MdlTypesNames;
MdlOrder = p.Results.MdlOrder;
SigLvl = p.Results.SigLvl;

pos = get(gca,'position');

ModelColors = MdlTypesCmap(MdlTypesInd,:);

mean_val = zeros(Ngrps,Nmodels);

% Find the order to plot the models in (if not preselected). Start by 
if isa(MdlOrder,'double')
    if length(MdlOrder) ~= Nmodels
        error('''order'' needs to specificy the ordering for all models')
    else
        mdl_order = MdlOrder;
    end
else
    
for i = 1:length(MdlNames)
    for j = 1:Ngrps
    mean_val(j,i) = mean(MdlData{j,i});
    
    if j == 2
    
        if isnan(mean_val(1,i)) && isnan(mean_val(2,i))
           error(['Model index number ',num2str(i),' has no data in either group!'])  
        end

        if isnan(mean_val(1,i))
            mean_val(1,i) = mean_val(2,i);
        end
        if isnan(mean_val(2,i))
            mean_val(2,1) = mean_val(1,i);
        end
        
    end
    
    end
end

switch MdlOrder
    case 'min'
    % Find which model type has the mode with the smallest mean value and
    % sort by that model type
        if Ngrps == 2
            [~,mdl_type_to_sort_by] = min(min(mean_val,[],2));
            [~,mdl_order] = sort(mean_val(mdl_type_to_sort_by,:),'ascend');
        elseif Ngrps == 1
            [~,mdl_order] = sort(mean_val,'ascend');    
        end
    case 'max'
    % Same as above but for the largest mean value
        if Ngrps == 2
            [~,mdl_type_to_sort_by] = max(max(mean_val,[],2));
            [~,mdl_order] = sort(mean_val(mdl_type_to_sort_by,:),'descend');
        elseif Ngrps == 1
            [~,mdl_order] = sort(mean_val,'descend');    
        end  
    otherwise
        error('Unknown order option specified')
end

end

ordered_mdl_data = MdlData(:,mdl_order);

%Set up VInd and xpos to iterate over to control the sequential ordering
%of figure handles to V and position of the violin plots respectively
VInd = 1;
xpos = 0;
for i = 1:Nmodels

    for j = 1:Ngrps
        xpos = xpos + 1;
        if ~isempty(ordered_mdl_data{j,i})
            % Checks if the corresponding subtype model has data
            % associated with it. If not, will centre the model which has
            % data when plotting
            if j == 1
                linestyle = '-';
                if Ngrps > 1
                    if isempty(ordered_mdl_data{2,i})
                        centre_violin = .5;
                    else
                        centre_violin = 0;
                    end
                else
                    centre_violin = 0;
                end
            elseif j == 2
                linestyle = ':';
                if isempty(ordered_mdl_data{1,i})
                    centre_violin = -.5;
                else
                    centre_violin = 0;
                end
            end

            V(1,VInd) = Violin(ordered_mdl_data{j,i}, xpos+centre_violin,'ShowMean',true,'ShowData',false,'ViolinAlpha',1);
            V(1,VInd).ViolinColor = ModelColors(mdl_order(i),:);
            V(1,VInd).MeanPlot.Color = [0 0 0];

            V(1,VInd).EdgeColor = [0 0 0];
            V(1,VInd).ViolinPlot.LineWidth = 2;
            V(1,VInd).ViolinPlot.LineStyle = linestyle;

            V(1,VInd).BoxColor = [0 0 0];

            VInd = VInd + 1;
        end

    end

end    

xlim([.5 (Nmodels*Ngrps)+.5])

% get the current ylimits as to reapply them later
ylimits = ylim;

xticklabels([]);

% Make new xtick positions for model subtype labels
if Ngrps == 2
    xticks(1.5:2:((Nmodels*2)+1));
elseif Ngrps == 1
    xticks(1:Nmodels)
end
xticklabels(MdlNames(mdl_order))
xtickangle(45)

ylabel(DataLabel)
xlabel('Models')
hold on
% As model subtypes are paired in each model, to make this more salient
% lines are plotted between different model pairs

if Ngrps == 2 

for i = 1:length(mdl_order)
    plot([-.5 -.5]+((2*mdl_order(i))-1),[-10 10],'Color','k','LineWidth',1.5);
    plot([.5 .5]+(2*mdl_order(i)),[-10 10],'Color','k','LineWidth',1.5);  
end

end

set(gca,'Fontsize',20)
set(gca, 'Layer','top'); set(gca,'TickDir','In')
ylim(ylimits)

if ~isempty(MdlTypesNames)

% To make a legend with the colour for each model type, get the
% number of types and then make a scatter plot just outside the FOV in the
% respective model type colour.
U = unique(MdlTypesInd);
for i = 1:length(U)
    lgd1_s(i) = scatter(2*length(mdl_order)+1,0,20,MdlTypesCmap(U(i),:),'filled');
    
end

% Make the legend for just the model types
lgd1 = legend(lgd1_s,MdlTypesNames,'Location','northoutside','Fontsize',14,'Orientation','horizontal','AutoUpdate','off');

end

if Ngrps == 2

% To make the legend for each model group, do the same as above but this
% time for plots with different line styles
lgd2_p(1) = plot([2*length(mdl_order)+1 2*length(mdl_order)+2],[1 1],'Color',[0 0 0],'LineWidth',2,'LineStyle','-');
lgd2_p(2) = plot([2*length(mdl_order)+1 2*length(mdl_order)+2],[1 1],'Color',[0 0 0],'LineWidth',2,'LineStyle',':');

%lgdPosition = lgd1.Position;

% Manually adjust the position of the legend
lgd1.Position(2) = 0.9171;

if ~isempty(MdlTypesNames)

% Plot significant stars over the models which have significant differents 
if SigLvl ~= 0
    if SigLvl < 1
        [issig,pval] = CompareFitsSigTest(ordered_mdl_data,SigLvl,0);
    elseif SigLvl == 1
        [issig,pval] = CompareFitsSigTest(ordered_mdl_data,.05,1);    
    elseif SigLvl > 1
        [issig,pval] = CompareFitsSigTest(ordered_mdl_data,.05,SigLvl);
    end
    
    for i = 1:Nmodels
        val_pairs{i} = [(i*2)-1 (i*2)];
    end
    
    sigstar(val_pairs(logical(issig)),pval(logical(issig)))
end

else
   issig = [];
   pval = [];
end

% Make a new invisible axis over the first. Needed because a plot cannot
% have two seperate legends
ah1=axes('position',pos,'visible','off');

% Make the legend for the model types
lgd2 = legend(ah1,lgd2_p,GrpNames,'Location','northoutside','Fontsize',14,'Orientation','horizontal','AutoUpdate','off');

%lgdPosition = lgd2.Position;

% Manually adjust its position
lgd2.Position(2) = 0.8622;

else
   issig = [];
   pval = [];    
end