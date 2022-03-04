function MakeFigure1A_and_2

% This makes Figure 1A and 2 for the paper

load('Fetal_faces_vertices_native_ss.mat')
load('fsaverage_surface_data.mat')

% Note the plotting of the figures took a bet of messing around to make it
% work appropriately and will probably break when tried on a different
% resolution screen (FYI this was done on a 2560*1440 screen). Basically I
% spent a lot of time making sure all the subplots were close enough
% together, then I use the exportgraphics option to trim all the white
% space
%
% I combined the images seperately in PowerPoint

figure('Position',[1 41  2560 1323])

lhsurface.faces = lh_faces;
lhsurface.vertices = lh_verts;
p = patch(lhsurface);
set(p,'EdgeColor','none','FaceColor','flat');
p.FaceLighting = 'gouraud';
lhsurface.vertices = fetal_vertices{18};
lhsurface.faces = fetal_faces{18}+1;
material dull
camlight;
camlight(-80,-10);
view([-90 0])
p = patch(lhsurface);
p.FaceColor = [.5 .5 .5];
% Manually set axis limits
xlimits = xlim;
ylimits = [-55 55];
zlimits = [-30 50];
clf

for i = 1:18
ax{i} = my_subplot(3,6,i,'OuterPosition');
lhsurface.vertices = fetal_vertices{i};
lhsurface.faces = fetal_faces{i}+1;

min_z = min(fetal_vertices{i}(:,3));

lhsurface.vertices(:,3) = lhsurface.vertices(:,3)+(-30-min_z);

p = patch(lhsurface);
set(p,'EdgeColor','none','FaceColor','flat');
p.FaceLighting = 'gouraud';
material dull
camlight;
camlight(-80,-10);
view([-90 0])
p.FaceColor = [.5 .5 .5];

axis off
axis tight
axis equal

xlim(xlimits)
ylim(ylimits)
zlim(zlimits)
ylabel(ax{i},['GA ',num2str(i+20)],'FontSize',32)

set(get(ax{i},'YLabel'),'Visible','on')

end

PosData = load('FetalPicPosition.mat');

for i = 1:18
    ax{i}.Position = PosData.Position(i,:);
end
exportgraphics(gcf,'FetalBrain.png','Resolution','300')

clf

for i = 1:18
ax{i} = my_subplot(3,6,i,'OuterPosition');
lhsurface.vertices = fetal_vertices{i};
lhsurface.faces = fetal_faces{i}+1;

min_z = min(fetal_vertices{i}(:,3));

lhsurface.vertices(:,3) = lhsurface.vertices(:,3)+(-30-min_z);

p = patch(lhsurface);
set(p,'FaceVertexCData',fetal_sulc{i},'EdgeColor','none','FaceColor','interp');
p.FaceLighting = 'gouraud';
material dull
camlight;
camlight(-80,-10);
view([-90 0])
colormap(turbo(256))
axis off
axis tight
axis equal

xlim(xlimits)
ylim(ylimits)
zlim(zlimits)
ylabel(ax{i},['GA ',num2str(i+20)],'FontSize',32)

set(get(ax{i},'YLabel'),'Visible','on')

end

for i = 1:18
    ax{i}.Position = PosData.Position(i,:);
end

exportgraphics(gcf,'FetalBrainSulc.png','Resolution','300')

load('fsaverage_surface_data.mat', 'lh_faces','lh_verts','lh_rand200')

% The surfaces we used above where in each fetal brains own native space.
% To plot the parcellation we need the surfaces aligned to fsaverage

load('Fetal_faces_vertices.mat')
for i = 1:18
ax{i} = my_subplot(3,6,i,'OuterPosition');
lhsurface.vertices = fetal_vertices_adult{i};
lhsurface.faces = fetal_faces_adult{i};

min_z = min(lhsurface.vertices(:,3));

lhsurface.vertices(:,3) = lhsurface.vertices(:,3)+(-30-min_z);

p = plotSurfaceROIBoundary(lhsurface,lh_rand200,1:100,'midpoint',lines(100),1,2);

p.FaceLighting = 'gouraud';
material dull
camlight;
camlight(-80,-10);
view([-90 0])
axis off
axis tight
axis equal

xlim(xlimits)
ylim(ylimits)
zlim(zlimits)
ylabel(ax{i},['GA ',num2str(i+20)],'FontSize',32)

set(get(ax{i},'YLabel'),'Visible','on')

end

for i = 1:18
    ax{i}.Position = PosData.Position(i,:);
end

exportgraphics(gcf,'FetalBrainParc.png','Resolution','300')
%% Make panels B and C

load('Fetal_faces_vertices_native_ss.mat')
load('fsaverage_surface_data.mat')

for i = 1:19

    if i < 19
lhsurface.vertices = fetal_vertices{i};
lhsurface.faces = fetal_faces{i}+1;
    else
lhsurface.vertices = lh_verts;
lhsurface.faces = lh_faces;        
    end

vert = lhsurface.vertices;
faces = lhsurface.faces;

a = vert(faces(:, 2), :) - vert(faces(:, 1), :);
b = vert(faces(:, 3), :) - vert(faces(:, 1), :);
c = cross(a, b, 2);
areaTri = 1/2 *(sqrt(sum(c.^2, 2)));

SA(i) = sum(areaTri);

end

figure('Position',[538 235 1440 639])

ax = axes('Position',[0.0681    0.1541    0.3966    0.7709]);
bar(SA)

ylabel('Surface area')

xlabel('Age (GA)');

% Onyl label certain weeks

GA_weeks = 1:3:16;

for i = 1:length(GA_weeks)
    cLabels{i} = num2str(GA_weeks(i)+20);    
end

cLabels{length(GA_weeks)+1} = 'Adult';

xticks(1:3:19)

xticklabels(cLabels);

set(gca, 'FontSize',18)
xtickangle(45)

ax = axes('Position',[0.5903    0.1541    0.3593    0.7709]);

load('random200_distances.mat')
%cmap = parula(19);

cmap = turbo(19);

for i = 1:19
distances = triu2vec(ADJS{i},1);

hold on

[f,x] = ksdensity(distances); 
plot(x,f,'Color', cmap(i,:),'LineWidth',2)
end

ylabel('Proportion of connections')
xlabel('Fibre distance (mm)')

colormap(cmap)

caxis([.5 19.5])

c = colorbar;

c.Label.String = 'Age (GA)';

GA_weeks = 1:3:16;

c.Ticks = 1:3:19;

for i = 1:6
cLabels{i} = [num2str(GA_weeks(i)+20)];

end
cLabels{7} = 'Adult';

c.TickLabels = cLabels;

set(gca, 'FontSize',20)

annotation('textbox',...
    [0.009629050341967,0.918918918918919,0.029776992936428,0.070945945945946],...
    'String',{'B'},...
    'LineStyle','none',...
    'FontSize',36,...
    'FitBoxToText','off');

annotation('textbox',...
    [0.496952461038233,0.918918918918919,0.029776992936428,0.070945945945946],...
    'String','C',...
    'LineStyle','none',...
    'FontSize',36,...
    'FitBoxToText','off');

saveas(gcf,'Figure2BC','svg')