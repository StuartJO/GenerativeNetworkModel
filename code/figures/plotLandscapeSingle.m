function plotLandscapeSingle(LANDSCAPE,X_points,Y_points,type,BestParams,ylabelname,xlabelname,xonly)

%labels = {'spatial','neighbors','matching','clu-avg','clu-min','clu-max','clu-diff','clu-prod','deg-avg','deg-min','deg-max','deg-diff','deg-prod','com'};



if nargin < 7
    xlabelname = '\eta';
end

if nargin < 8
    xonly = 0;
end

if xonly ~= 2
X_points = X_points*-1;
end
%x_points = x_points*-1;

lines_cmap = lines(4);
alphavals = linspace(0,1,256);
switch type 
        case 'deg'
            for i = 1:256
                cmap(i,:) = findAlphaColor(lines_cmap(1,:),alphavals(i));
            end
        case 'clu'
            for i = 1:256
                cmap(i,:) = findAlphaColor(lines_cmap(2,:),alphavals(i));
            end            
        case 'bet'
            for i = 1:256
                cmap(i,:) = findAlphaColor(lines_cmap(3,:),alphavals(i));
            end            
        case 'len'
            for i = 1:256
                cmap(i,:) = findAlphaColor(lines_cmap(4,:),alphavals(i));
            end            
end

    
    if xonly == 1
            fakeMat = [ones(1,length(X_points)); ones(1,length(X_points))*2; ones(1,length(X_points))*3];
            
                switch type
                    case 'lim'
                        meanLand = mode(LANDSCAPE);
                    otherwise
                        meanLand = nanmean(LANDSCAPE);
                end
            
            
            s = pcolor(X_points(1:3,:),fakeMat,[meanLand; meanLand; meanLand]);  
            yticks([])
            
    elseif xonly == 2
            fakeMat = [ones(length(Y_points),1) ones(length(Y_points),1)*2 ones(length(Y_points),1)*3];
            
                switch type
                    case 'lim'
                        meanLand = mode(LANDSCAPE,2);
                    otherwise
                        meanLand = nanmean(LANDSCAPE,2);
                end
            
            
            s = pcolor(fakeMat,Y_points(:,1:3),[meanLand meanLand meanLand]);  
            xticks([])            
    else
    
    s = pcolor(X_points,Y_points,LANDSCAPE);
    
    end
        
    s.LineStyle = 'none';
    
    switch type
        case 'corr'
            s.FaceColor = 'interp';
            caxis([-.5 .5])
            colormap([make_cmap('steelblue',50,30,0);flipud(make_cmap('orangered',50,30,0))])
        case 'ks'
            s.FaceColor = 'interp';
            caxis([0 1])
        case 'lim'
            colormap(lines(4))
            caxis([1 4])
        case {'deg','bet','clu','len'}
            s.FaceColor = 'interp';
            caxis([0 1])            
            colormap(cmap)      
    end
    hold on
    
    if ~isempty(BestParams)
    if xonly == 1
    scatter(BestParams(:,1)*-1,(BestParams(:,2).^0)*2,20,'k+')   
    elseif xonly == 2
    scatter((BestParams(:,2).^0)*2,BestParams(:,1),20,'k+')       
    else
    scatter(BestParams(:,1)*-1,BestParams(:,2),20,'k+')
    end
    end
    
ylabel(ylabelname)
xlabel(xlabelname)

ax = gca;

set(ax,'Fontsize',16)
ax.TickLength = [.025, .025]; % Make tick marks longer.
ax.LineWidth = 100*.01; % Make tick marks thicker.
set(gca, 'Layer', 'top')

% 
% c = colorbar('Position',[0.9212    0.0794    0.0171    0.8450]);
% 
% switch type
%     case 'corr'
%        c.Label.String = 'Spearman correlation with emperical degree';
%     case 'ks'
%        c.Label.String = 'Spearman correlation with emperical degree';
% end
% 
% c.Label.String = 'Spearman correlation with emperical degree';
% c.LineWidth = 2;
% set(gca,'Fontsize',28)
run = 0;
if run == 1
    switch type
        case 'corr'
colormap([make_cmap('steelblue',256,30,0);flipud(make_cmap('orangered',256,30,0))])
c = colorbar('Position',[0.1293    0.068    0.7758    0.035],'Orientation','Horizontal');
c.Label.String = 'Degree rank correlation';
c.FontSize = 20;
c.Label.FontSize = 28;
caxis([-.5 .5])
c.LineWidth = 2;
%set(c, 'xlim', [0 1])
%colormap(flipud(make_cmap('orangered',100,30,0)))

%axis off

        case 'ks'
            c = colorbar('Position',[0.1293    0.068    0.7758    0.035],'Orientation','Horizontal');
            c.Label.String = 'max(\it{KS})';
            c.FontSize = 20;
            c.Label.FontSize = 28;
            caxis([0 1])
            set(c, 'xlim', [0 1])
            c.LineWidth = 2;
            %colormap(flipud(make_cmap('orangered',100,30,0)))


    
    
    case 'lim'
            colormap(lines(4))
    c = colorbar('Position',[0.1293    0.068    0.7758    0.035],'Orientation','Horizontal','Ticks',1:4,'TickLabels',{'Degree','Clustering','Betweenness','Edge Length'});
    %c.Label.String = 'max(\it{KS})';
    c.FontSize = 20;
    c.Label.FontSize = 28;
    caxis([0.5 4.5])
    c.LineWidth = 2;
    %set(c, 'xlim', [0 1])
    %colormap(flipud(make_cmap('orangered',100,30,0)))
    end
    
end