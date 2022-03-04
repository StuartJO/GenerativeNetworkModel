function MakeFibreDistanceExampleFig

% This will plot Figure S7 from the paper. It will actually perform a 2D
% version of what is proposed in the paper

%Equation for semi circle
    th = linspace( pi, 0, 99);
    R = 1;  %or whatever radius you want
    x_out = R*cos(th);
    y_out = R*sin(th);
    y_out(1) = 0;

%Equation for semi circle
    th = linspace( 0, pi, 59);
    R = .7;  %or whatever radius you want
    x_in = R*cos(th);
    y_in = R*sin(th);
    y_in(1) = 0;

    n = 7;
    
    right_bottom_x = linspace(1,.7,n); 
    right_bottom_y = zeros(1,n);
    
    left_bottom_x = fliplr(right_bottom_x*-1);
    left_bottom_y = right_bottom_y;
    
    
    allx_line = [x_out right_bottom_x(2:n-1) x_in left_bottom_x(2:n)];
    ally_line = [y_out right_bottom_y(2:n-1) y_in left_bottom_y(2:n)];
    
%     plot(allx_line,ally_line)    
    
    redu_out = round(linspace(1,99,15));
    redu_in = round(linspace(1,59,9));
    
    redu_x_out = x_out(redu_out);
    redu_y_out = y_out(redu_out);
    
    redu_x_in = x_in(redu_in);
    redu_y_in = y_in(redu_in);
    

    
    redu_x_line = [x_out(redu_out) x_in(redu_in) x_out(1)];
    redu_y_line = [y_out(redu_out) y_in(redu_in) y_out(1)]; 
%     hold on
%     plot(redu_x_line,redu_y_line)
        
%     plot(x_in(redu_in)*1.05,y_in(redu_in)+.015)
%     
%     plot(x_out(redu_out)*.95,y_out(redu_out)-.01)
    
    for i = 1:length(x_out)
        if i == 1
        mid_point_x_out(i) = (left_bottom_x(n-1)+x_out(i+1))/2;
        mid_point_y_out(i) = (left_bottom_y(n-1)+y_out(i+1))/2;            
        elseif i == length(x_out)
        mid_point_x_out(i) = (x_out(i-1)+right_bottom_x(2))/2;
        mid_point_y_out(i) = (y_out(i-1)+right_bottom_y(2))/2;     
        else
        mid_point_x_out(i) = (x_out(i-1)+x_out(i+1))/2;
        mid_point_y_out(i) = (y_out(i-1)+y_out(i+1))/2;
        end
        
    end

    for i = 1:length(x_in)
        if i == 1
        mid_point_x_in(i) = (right_bottom_x(n-1)+x_in(i+1))/2;
        mid_point_y_in(i) = (right_bottom_y(n-1)+y_in(i+1))/2;            
        elseif i == length(x_in)
        mid_point_x_in(i) = (x_in(i-1)+left_bottom_x(2))/2;
        mid_point_y_in(i) = (y_in(i-1)+left_bottom_y(2))/2;     
        else
        mid_point_x_in(i) = (x_in(i-1)+x_in(i+1))/2;
        mid_point_y_in(i) = (y_in(i-1)+y_in(i+1))/2;
        end
        
    end
    
    mid_point_redu_out= [mid_point_x_out(redu_out)' mid_point_y_out(redu_out)'];
    mid_point_redu_in= [mid_point_x_in(redu_in)' mid_point_y_in(redu_in)'];
    
    Incoords = find_point_on_line([redu_x_in' redu_y_in'],mid_point_redu_in,-.02,'distance');
    
    Incoords(1,:) = find_point_on_line([redu_x_in(1) redu_y_in(1)],mid_point_redu_in(1,:),.02,'distance');
    Incoords(end,:) = find_point_on_line([redu_x_in(end) redu_y_in(end)],mid_point_redu_in(end,:),.02,'distance');
    
    Outcoords = find_point_on_line([redu_x_out' redu_y_out'],mid_point_redu_out,.02,'distance');
    
    Subsurf_points = [[Outcoords(:,1); Incoords(:,1)],[Outcoords(:,2); Incoords(:,2)]];
    
    [distance, direct_connection, eucl_dist] = distance_through_shape(Subsurf_points(:,1),Subsurf_points(:,2),0,allx_line,ally_line,1);
        
    G = graph(direct_connection.*distance);
    
    SubSurfColor = [247 46 70]./255;  

    fillColor = [199 212 237]./255;
    figure('Position',[489 367 2018 849])
    tiledlayout(2,3,'TileSpacing','none','Padding','none')
    nexttile
    
    fill(allx_line,ally_line,fillColor)
    hold on
    scatter(allx_line(1:end-1),ally_line(1:end-1),'k','filled')
    axis off
    Lbl = AddPlotLabel(gca,'A',40,0,.425,2);
    nexttile
    fill(redu_x_line,redu_y_line,fillColor)
    hold on
    scatter(redu_x_line(1:end-1),redu_y_line(1:end-1),'k','filled')
    axis off
    Lbl = AddPlotLabel(gca,'B',40,0,.425,2);
    nexttile
    fill(redu_x_line,redu_y_line,fillColor)
    hold on
    scatter(redu_x_line(1:end-1),redu_y_line(1:end-1),'k','filled')
    scatter(Subsurf_points(:,1),Subsurf_points(:,2),'MarkerEdgeColor','k','MarkerFaceColor',SubSurfColor)
    axis off
    Lbl = AddPlotLabel(gca,'C',40,0,.425,2);
    nexttile
    fill(allx_line,ally_line,[199 212 237]./255)
    hold on
    %scatter(redu_x_line(1:end-1),redu_y_line(1:end-1),'k','filled')
    scatter(allx_line(1:end-1),ally_line(1:end-1),'k','filled')
    axis off
    
    plot([Subsurf_points(1,1),Subsurf_points(7,1)],[Subsurf_points(1,2),Subsurf_points(7,2)],'LineWidth',4,'Color',[146 208 80]./255)
    plot([Subsurf_points(9,1),Subsurf_points(11,1)],[Subsurf_points(9,2),Subsurf_points(11,2)],'LineWidth',4,'Color',[146 208 80]./255)
    plot([Subsurf_points(13,1),Subsurf_points(17,1)],[Subsurf_points(13,2),Subsurf_points(17,2)],'LineWidth',4,'Color',[146 208 80]./255)
    
    plot([Subsurf_points(23,1),Subsurf_points(12,1)],[Subsurf_points(23,2),Subsurf_points(12,2)],'LineWidth',4,'Color','r')
    plot([Subsurf_points(21,1),Subsurf_points(19,1)],[Subsurf_points(21,2),Subsurf_points(19,2)],'LineWidth',4,'Color','r')
    plot([Subsurf_points(24,1),Subsurf_points(16,1)],[Subsurf_points(24,2),Subsurf_points(16,2)],'LineWidth',4,'Color','r')
        scatter(Subsurf_points(:,1),Subsurf_points(:,2),'MarkerEdgeColor','k','MarkerFaceColor',SubSurfColor)
    Lbl = AddPlotLabel(gca,'D',40,0,.425,2);
    nexttile
    f = plot(G);
    f.NodeLabel = [];
    %     
     f.XData = redu_x_line(1:end-1);
     f.YData = redu_y_line(1:end-1);
     f.NodeLabel = [];
     f.LineStyle = ':';
     f.EdgeColor = 'k';
     f.NodeColor = 'k';
     f.LineWidth = 1;
     xlim([-1 1])
     ylim([0 1])
     axis off
     hold on
     scatter(redu_x_line(1:end-1),redu_y_line(1:end-1),'k','filled')
     Lbl = AddPlotLabel(gca,'E',40,0,.425,2);
     
     nexttile
    fill(redu_x_line,redu_y_line,[199 212 237]./255)
    hold on
    %scatter(redu_x_line(1:end-1),redu_y_line(1:end-1),'k','filled')
    axis off
    
        f = plot(G);
    f.NodeLabel = [];
    %     
     f.XData = redu_x_line(1:end-1);
     f.YData = redu_y_line(1:end-1);
     f.NodeLabel = [];
     f.EdgeAlpha = 0;
     f.EdgeColor = 'k';
     f.NodeColor = 'k';
     xlim([-1 1])
     ylim([0 1])
     axis off
    
    scatter(redu_x_line(1:end-1),redu_y_line(1:end-1),'k','filled')
    path = shortestpath(G,1,13);
    plot(redu_x_line(path),redu_y_line(path),':','LineWidth',4,'Color',[146 208 80]./255)
    plot(redu_x_line([1 13]),redu_y_line([1 13]),':','LineWidth',4,'Color','r')
    Lbl = AddPlotLabel(gca,'F',40,0,.425,2);
    
    set(0,'defaultAxesFontName', 'Arial')
    set(0,'defaultTextFontName', 'Arial')
    
    %SaveFigures('Chap5_Fig2_take2',1);
    %saveas(gcf,'Chap5_Fig2','svg')
    %saveas(gcf,'FibreDistanceExample','svg')