function PlotModSurf(Emp,Model,Model_name,plotLabel)

Emp = [Emp(1:53) NaN Emp(54:99)];

Model = [Model(1:53) NaN Model(54:99)];
clim_max = max(max(Model),max(Emp));
climits1 = [1 clim_max];
climits2 = [1 clim_max];

figure('Position',[0 0 1035 917])

load('fsaverage_surface_data.mat','lh_inflated_verts','lh_faces','lh_rand200')
%addpath ./cbrewer

%lh_rand200(lh_rand200==54) = 0;

% cmap = flipud(make_cmap('orangered',100,30,0));
cmap = lines(clim_max);

lhsurface.vertices = lh_inflated_verts;
lhsurface.faces = lh_faces;

latpatch1 = axes('Position',[0.0100    0.5650    0.4315    0.4762]);
plotSurfaceROIBoundary(lhsurface,lh_rand200,Emp,'midpoint',cmap,1,2,climits1);

camlight(80,-10);
camlight(-80,-10);
camlight(180,0);
view([-90 0])
axis off
axis tight
axis equal

medpatch1 = axes('Position',[0.4561    0.5650    0.4315    0.4762]);
plotSurfaceROIBoundary(lhsurface,lh_rand200,Emp,'midpoint',cmap,1,2,climits1);
camlight(80,-10);
camlight(-80,-10);
camlight(180,0);
view([90 0])
axis off
axis tight
axis equal

latpatch2 = axes('Position',[0.0100    0.1450    0.4315    0.4762]);
plotSurfaceROIBoundary(lhsurface,lh_rand200,Model,'midpoint',cmap,1,2,climits2);

camlight(80,-10);
camlight(-80,-10);
camlight(180,0);
view([-90 0])
axis off
axis tight
axis equal


medpatch2 = axes('Position',[0.4561    0.1450    0.4315    0.4762]);
plotSurfaceROIBoundary(lhsurface,lh_rand200,Model,'midpoint',cmap,1,2,climits2);
camlight(80,-10);
camlight(-80,-10);
camlight(180,0);
view([90 0])
axis off
axis tight
axis equal
%0.8776
TEXTBOX_TOP = annotation('textbox',[0    0.5530 0.9077 0.1093],'String','Empirical','EdgeColor','none','FontSize',36,'HorizontalAlignment','center');

TEXTBOX_BOTTOM = annotation('textbox',[0    0.1124  0.9077 0.1093],'String',Model_name,'EdgeColor','none','FontSize',36,'HorizontalAlignment','center');

%ColorbarLim = climits;

% c1 = colorbar(medpatch1,'Location','West');
% c1.AxisLocation = 'out';
% c1.Position = [0.8947 0.6472 0.0251 0.3];
% c1.Label.String = 'Degree';
% c1.FontSize = 20;
% c1.TickLength = 0.04;
% %c1.Label.Position = [2.7619   15.9900 0];
% set(c1, 'xlim', climits1)
% 
% 
% c2 = colorbar(medpatch2,'Location','West');
% c2.AxisLocation = 'out';
% c2.Position = [0.8947 0.2072 0.0251 0.3];
% c2.Label.String = 'Degree';
% c2.FontSize = 20;
% c2.TickLength = 0.04;
% %c2.Label.Position = [2.7619   15.9900 0];
% set(c2, 'xlim', climits2)


annotation('textbox',[0.0117    0.9039    0.0588    0.0885],'String',plotLabel,'EdgeColor','none','FontSize',48);
