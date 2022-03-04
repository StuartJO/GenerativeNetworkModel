function VisualiseKSstats(KS,MdlLabels,NPlotrows,NPlotcolumns)

% This function will plot the distribution of KS values for different 
% models, and label the proportion of times each of the four properties
% determines maxKS
%
% Input:
%
% KS = an 1*nMdls cell, each containing an N*4 matrix. Each row corresponds
% to the KS values for a particular network, each column corresponds to a 
% particular KS measure in the following order: 'Degree','Clustering',
% 'Betweenness','Edge length'
%
% MdlLabels = a cell containing the name of each model
%
% NPlotrows = the number of rows in the plot
%
% NPlotcolumns = the number of columns in the plot 
% 

PlotColors = [0.1103    0.4122    0.6888;...
    0.8099    0.2586         0;...
    0.1086    0.5167    0.2470;...
    0.3949    0.2767    0.6193];

figure('Position',[1 41 1920 963])
t = tiledlayout(NPlotrows,NPlotcolumns);
for i = 1:length(KS)
       
    maxKS = zeros(size(KS{i},1),4);
    
    % Indicate with a binary value which of the 4 properties determined
    % maxKS
    
    [~,maxI] = max(KS{i},[],2,'linear');
    maxKS(maxI) = 1;
    
    % Take the mean, find the proportion
    
    maxKS_data = mean(maxKS);
    
   nexttile
   
   % Plot the distributions
   
   for j = 1:4
        V = Violin(KS{i}(:,j), j,'ShowMean',true,'ShowData',false,'ViolinAlpha',1);
        V.ViolinColor = PlotColors(j,:);
        V.MeanPlot.Color = [0 0 0];
        V.EdgeColor = [0 0 0];
        V.ViolinPlot.LineWidth = 1;
        V.BoxColor = [0 0 0];
   end
   
    % Add titles
    title(MdlLabels{i})
    ylimits = ylim;
    yticks('manual')
    xlim([.5 4.5])
    % Extend the ylimits a small amount so the proportions can be written
    % above each distribution
    y_max = find_point_on_line(ylimits(1),ylimits(2),1.12);
    ylim([ylimits(1) y_max])
    % Find the y position to write the proportions
    y_text = find_point_on_line(ylimits(1),ylimits(2),1.07);
    % Write the proportions above each plot
    for j = 1:4
        text(j,y_text,[num2str(round(maxKS_data(j)*100,1)),'%'],'HorizontalAlignment','center','FontSize',14)
    end
    set(gca,'FontSize',18)
    xticks([])
end

% Make some points outside of the axis limits. Use these to make the legend
for i = 1:4
    s(i) = scatter(5,0,20,PlotColors(i,:),'filled');
    hold on
end

% Make a fake axis to assign the legend to so it is easier to place
pos = get(gca,'position');
ah1=axes('position',pos,'visible','off');

lgd1 = legend(ah1,s,{'Degree','Clustering','Betweenness','Edge length'},'Location','northoutside','Fontsize',20,'Orientation','horizontal','AutoUpdate','off');

% Manually reposition the legend 
lgd1.Position(2) = 0.002;
lgd1.Position(1) = (1-(lgd1.Position(3)))/2;
t.TileSpacing = 'compact';
t.Padding = 'compact';
ylabel(t,'KS statistic','FontSize',28)