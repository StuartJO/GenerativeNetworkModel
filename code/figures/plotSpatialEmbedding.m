function plotSpatialEmbedding(Emp,Model,D,MissingROI)

% This function plots the CDF

figure('Position',[1 41 1920 963])

if iscell(Emp)

    for i = 1:length(Emp)
        
        A = double(Emp{i}>0);
        
        EmpData{1}(i,:) = sum(A); 


        EmpData{3}(i,:) = betweenness_bin(A);


        EmpData{2}(i,:) = clustering_coef_bu(A)';

        EmpData{4}(i,:) = sum(D.*A)/2; 

    end
    
        for i = 1:length(Model)
        ModelData{1}(i,:) = sum(Model{i});    


        ModelData{3}(i,:) = betweenness_bin(Model{i});  


        ModelData{2}(i,:) = clustering_coef_bu(Model{i})';    

        ModelData{4}(i,:) = sum(D.*Model{i})/2;    

        end

else
    i = 1;
    Emp = double(Emp>0);
        EmpData{1}(i,:) = sum(Emp);
        ModelData{1}(i,:) = sum(Model);    


        EmpData{3}(i,:) = betweenness_bin(Emp);
        ModelData{3}(i,:) = betweenness_bin(Model);  


        EmpData{2}(i,:) = clustering_coef_bu(Emp)';
        ModelData{2}(i,:) = clustering_coef_bu(Model)';    

        EmpData{4}(i,:) = sum(D.*Emp)/2;
        ModelData{4}(i,:) = sum(D.*Model)/2;       
end

lines_cmap = lines(4);

for i = 1:4
%lines_cmap_alpha{i} = make_alpha_cmap(lines_cmap(i,:),256);
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
lines_cmap_alpha{i} = make_alpha_cmap(lines_cmap_alpha{i}(256,:),256);
end

cbar_title = {'Degree','Clustering','Betweenness','Mean Edge length'};


for i = 1:4
    
    if i == 1
        row_offset = 0;
    elseif i == 2
        row_offset = 3;
    elseif i == 3
        row_offset = 12;
    elseif i == 4
        row_offset = 15;
    end
    
        if i == 1

    titlename = 'Degree';

elseif i == 2
    titlename = 'Clustering';
    
elseif i == 3
 
    
    titlename = 'Betweenness';

elseif i == 4

    titlename = 'Mean distance';
    end

CDFplotInd = 1+row_offset;
scatterplotInd = 7+row_offset;
medial_emp = 2+row_offset;
lateral_emp = 3+row_offset; 
medial_mdl = 8+row_offset;
lateral_mdl = 9+row_offset; 

my_subplot(4,6,CDFplotInd);

plotCDFsSingle(Emp,Model,D,i)
title(titlename)
Data1 = mean(EmpData{i},1);
Data2 = mean(ModelData{i},1);
my_subplot(4,6,scatterplotInd);
    scatter(Data1,Data2,20,lines_cmap_alpha{i}(128,:),'filled')
    
    xlabel('Observed')
    ylabel('Simulated')
    

    [c,pVal] = corr(Data1',Data2');

    hold on
    
    p = polyfit(Data1,Data2,1);
    f = polyval(p,Data1);
    plot(Data1,f,'Color',scatter_cmap(i,:),'LineWidth',3)
    
    
    xlimits = xlim;
    ylimits = ylim;
    
    text_loc_x = find_point_on_line(xlimits(1),xlimits(2),.65);
    text_loc_y = find_point_on_line(ylimits(1),ylimits(2),.8);
       
    text(text_loc_x,text_loc_y,{['{\itr} = ',num2str(round(c,3))],['{\itp} = ',num2str(round(pVal,3))]},'FontSize',12)
    set(gca,'FontSize',12)
    
    if ~isempty(MissingROI)
        Data1_allrois = [Data1(1:MissingROI-1) NaN Data1(MissingROI:end)];
        Data2_allrois = [Data2(1:MissingROI-1) NaN Data2(MissingROI:end)];
    else
        Data1_allrois = Data1;
        Data2_allrois = Data2;
    end
    
    climits = [0 nanmax(nanmax(Data1_allrois),nanmax(Data2_allrois))];
    ax1 = my_subplot(4,6,medial_emp);
    ax1.Position(1) = ax1.Position(1)-(0.0344/2);
    PlotPropSurfSingle(Data1_allrois,1,climits,lines_cmap_alpha{i})
    ax2 = my_subplot(4,6,lateral_emp);
    ax2.Position(1) = ax2.Position(1)-(0.0344);
    PlotPropSurfSingle(Data1_allrois,2,climits,lines_cmap_alpha{i})
    ax3 = my_subplot(4,6,medial_mdl);
    ax3.Position(1) = ax3.Position(1)-(0.0344/2);
    PlotPropSurfSingle(Data2_allrois,1,climits,lines_cmap_alpha{i})
    ax4 = my_subplot(4,6,lateral_mdl);
    ax4.Position(1) = ax4.Position(1)-(0.0344);
    PlotPropSurfSingle(Data2_allrois,2,climits,lines_cmap_alpha{i})
    
    cbar_x_start = ax4.Position(1)+ax4.Position(3)+.005;
    cbar_x_length = .01;
    cbar_y_start = ax4.Position(2);  
    cbar_y_length = [ax2.Position(2)+ax2.Position(4)]-ax4.Position(2);
    
    TEXTBOX_TOP = annotation('textbox',[ax1.Position(1) ax1.Position(2) ax2.Position(1)+ax2.Position(3)-ax1.Position(1) 0.05],'String','Empirical','EdgeColor','none','FontSize',36,'HorizontalAlignment','center');

    TEXTBOX_BOTTOM = annotation('textbox',[ax3.Position(1) ax3.Position(2) ax4.Position(1)+ax4.Position(3)-ax3.Position(1) 0.05],'String','Model','EdgeColor','none','FontSize',36,'HorizontalAlignment','center');

    
    c = colorbar;
    
    c.Position = [cbar_x_start cbar_y_start cbar_x_length cbar_y_length];
    c.Label.String = titlename;
    c.FontSize = 14;
    set(c, 'ylim', climits)
    
end

PanelA = annotation('textbox',[0 .97 .025 .025],'String','A','EdgeColor','none','FontSize',36,'HorizontalAlignment','center','VerticalAlignment','middle');
PanelB = annotation('textbox',[0.5 .97 .025 .025],'String','B','EdgeColor','none','FontSize',36,'HorizontalAlignment','center','VerticalAlignment','middle');
PanelC = annotation('textbox',[0 .47 .025 .025],'String','C','EdgeColor','none','FontSize',36,'HorizontalAlignment','center','VerticalAlignment','middle');
PanelD = annotation('textbox',[0.5 .47 .025 .025],'String','D','EdgeColor','none','FontSize',36,'HorizontalAlignment','center','VerticalAlignment','middle');
