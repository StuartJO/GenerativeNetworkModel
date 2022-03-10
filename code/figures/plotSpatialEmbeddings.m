function plotSpatialEmbeddings(EmpData,ModelData,CorrDATA,SurfaceData,MdlNames,MdlColors,SAVELOC)

% This function plots the distribution of correlations for a number of models,
% and the CDF and scatter plots and spatial projection of one specific model
%
% Note this script uses a bit of manual repositioning of axes

Nmodels = length(ModelData);

lines_cmap = lines(4);

for i = 1:4
if i == 1
lines_cmap_alpha{i} = brewermap(256,'Blues');
elseif i == 2
lines_cmap_alpha{i} = brewermap(256,'Oranges');     
elseif i == 3
lines_cmap_alpha{i} = brewermap(256,'Greens');     
elseif i == 4
lines_cmap_alpha{i} = brewermap(256,'Purples');    
end
scatter_cmap(i,:) = lines_cmap_alpha{i}(256,:);
line_cmap(i,:) = lines_cmap_alpha{i}(200,:);
end

cbar_title = {'Degree','Clustering','Betweenness','Mean Edge length'};

for i = 1:4
    for n = 1:length(ModelData)
            DATA{i,n} = CorrDATA{n}(:,i);
    end
end

for i = 1:4
        figure('Position',[1 41 1920 519])
        if i == 1

    titlename = 'Degree';
    lowercasename = 'degree';

elseif i == 2
    titlename = 'Clustering';
    lowercasename = 'clustering';
    
elseif i == 3
 
    
    titlename = 'Betweenness';
    lowercasename = 'betweenness';

elseif i == 4

    titlename = 'Mean distance';
    lowercasename = 'mean distance';
        end
    
 
my_subplot(1,2,1);

    PLOTIND = 1;
       for j = 1:size(DATA,2)
    V(1,PLOTIND) = Violin(DATA{i,j}, PLOTIND,'ShowMean',true,'ShowData',false,'ViolinAlpha',1);
    V(1,PLOTIND).ViolinColor = MdlColors(j,:);
    V(1,PLOTIND).MeanPlot.Color = [0 0 0];
    V(1,PLOTIND).BoxColor = [0 0 0];
    
    XtickLabels{PLOTIND} = MdlNames{j}; 
    PLOTIND = PLOTIND + 1;
       end    
     
xlim([.5 Nmodels+.5])

ylimits = ylim;

xticklabels([]);
xticks(1:Nmodels)

xticklabels(XtickLabels)

xtickangle(45)

ylabel({'Spearman correlation',['with ',lowercasename]})

xlabel('Model')

ylim(ylimits)

set(gca,'Fontsize',20)
set(gca, 'Layer','top'); set(gca,'TickDir','In')
ylim(ylimits)
        

CDFplotInd = 4;
scatterplotInd = 10;
medial_emp = 5;
lateral_emp = 6; 
medial_mdl = 11;
lateral_mdl = 12; 

cdfplot = my_subplot(2,6,CDFplotInd);
offset = .03;
plotCDFs(EmpData,ModelData,i,line_cmap)
title(titlename)

cdfplot.Position(1) = cdfplot.Position(1)-offset-.0119;
cdfplot.Position(3) = cdfplot.Position(3)+.0119;
cdfplot.Position(2) = cdfplot.Position(2)+.02;

Data1 = nanmean(EmpData{i},1);
X = ModelData{i};
X(isnan(X)) = 0;
Data2 = nanmean(X,1);

scatterplot = my_subplot(2,6,scatterplotInd);
scatterplot.Position(1) = scatterplot.Position(1)-offset;
scatterplot.Position(4) = scatterplot.Position(4)-.05;
    scatter(Data1,Data2,20,lines_cmap_alpha{i}(128,:),'filled')
    
    xlabel(['Empirical ',lowercasename])
    ylabel(['Model ',lowercasename])
    
    
    scatterplot.Position(2) = .125;
    %scatterplot.Position(4) = scatterplot.Position(4)-.119;
    
    [c,pVal] = corr(Data1',Data2','Type','Spearman');

    hold on
    
    % Fit a linear trend
    p = polyfit(Data1,Data2,1);
    f = polyval(p,Data1);
    plot(Data1,f,'Color',scatter_cmap(i,:),'LineWidth',3)
    
    
    xlimits = xlim;
    ylimits = ylim;
    ylim([ylimits(1) find_point_on_line(ylimits(1),ylimits(2),1.15)]);
    ylimits = ylim;
    text_loc_x = find_point_on_line(xlimits(1),xlimits(2),.5);
    text_loc_y = find_point_on_line(ylimits(1),ylimits(2),.95);
     
    set(gca,'FontSize',16)
    text(text_loc_x,text_loc_y,['\fontsize{18}{\rho}\fontsize{14} = ',num2str(round(c,3)),', {\itp} = ',num2str(round(pVal,3))],'HorizontalAlignment','center')
    
    xlim(xlimits)
    ylim(ylimits)
    
    if ~isempty(SurfaceData.MissingROI)
        Nnodes = max(SurfaceData.parc);
        MissingROI_ind = true(Nnodes,1);
        MissingROI_ind(SurfaceData.MissingROI) = false;
        Data1_allrois = nan(1,Nnodes);
        Data2_allrois = nan(1,Nnodes);
        Data1_allrois(MissingROI_ind) = Data1;
        Data2_allrois(MissingROI_ind) = Data2;
    else
        Data1_allrois = Data1;
        Data2_allrois = Data2;
    end
    climits = [0 nanmax(nanmax(Data1_allrois),max(nanmax(Data2_allrois)))];
    ax1 = my_subplot(2,6,medial_emp);
    ax1.Position(1) = ax1.Position(1)-(0.0344/2)-offset;
    PlotSurfaceData(Data1_allrois,1,SurfaceData,climits,lines_cmap_alpha{i})
    ax2 = my_subplot(2,6,lateral_emp);
    ax2.Position(1) = ax2.Position(1)-(0.0344)-offset;
    PlotSurfaceData(Data1_allrois,2,SurfaceData,climits,lines_cmap_alpha{i})
    ax3 = my_subplot(2,6,medial_mdl);
    ax3.Position(1) = ax3.Position(1)-(0.0344/2)-offset;
    ax3.Position(2) = 0.0858;
    PlotSurfaceData(Data2_allrois,1,SurfaceData,climits,lines_cmap_alpha{i})
    ax4 = my_subplot(2,6,lateral_mdl);
    ax4.Position(1) = ax4.Position(1)-(0.0344)-offset;
    ax4.Position(2) = 0.0858;
    PlotSurfaceData(Data2_allrois,2,SurfaceData,climits,lines_cmap_alpha{i})
    
    cbar_x_start = ax4.Position(1)+ax4.Position(3)+.005;
    cbar_x_length = .02;
    cbar_y_start = ax4.Position(2);  
    cbar_y_length = [ax2.Position(2)+ax2.Position(4)]-ax4.Position(2);
    
    TEXTBOX_TOP = annotation('textbox',[ax1.Position(1) ax1.Position(2) ax2.Position(1)+ax2.Position(3)-ax1.Position(1) 0.05],'String','Empirical','EdgeColor','none','FontSize',36,'HorizontalAlignment','center');

    TEXTBOX_BOTTOM = annotation('textbox',[ax3.Position(1) ax3.Position(2) ax4.Position(1)+ax4.Position(3)-ax3.Position(1) 0.05],'String','Model','EdgeColor','none','FontSize',36,'HorizontalAlignment','center');

    
    c = colorbar;
    
    c.Position = [cbar_x_start cbar_y_start cbar_x_length cbar_y_length];
    c.Label.String = titlename;
    c.FontSize = 18;
    set(c, 'ylim', climits)
    
    PanlLabls = {'A','B','C','D'};
    PanelLabel = annotation('textbox',[0 .95 .025 .025],'String',PanlLabls{i},'EdgeColor','none','FontSize',36,'HorizontalAlignment','center','VerticalAlignment','middle');   
    print([SAVELOC,PanlLabls{i},'.png'],'-dpng','-r300')
end