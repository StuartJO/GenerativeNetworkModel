function [p,h,c1,c2] = plotBrainProjectedConnectome(A,coords,NodeValues,options)

% This function makes a plot where the connectome is projected into 3D
% space with an overlay of the brain to aid with interpretation

% The outputs allow you to manipulate the figure after it has been made 
% WITH THE EXCEPTION OF EDGE AND NODE COLORS. Because of the way color 
% mapping works in this script, nodes/edges are applied an RGB value
% directly from a colormap, instead of node/edge CData being mapped onto a
% global colormap. However remaking the figure shouldn't be too hard as it
% renders pretty quick!

if nargin < 3
    NodeValues = [];
end

if nargin < 4
    options = struct;
end

if ~isfield(options,'NodeSize')
    options.NodeSize = 10;
end

if ~isfield(options,'EdgeColor')
    options.EdgeColor = [];
end 

if ~isfield(options,'NodeColor')
    if isempty(NodeValues)
        options.NodeColor = [0 0.4470 0.7410];
    else        
        options.NodeColor = [];
    end
end 

if ~isfield(options,'EdgeAlpha')
    options.EdgeAlpha = 1;
end 

if ~isfield(options,'EdgeSize')
    options.EdgeSize = 5;
end 

if ~isfield(options,'SurfaceAlpha')
    options.SurfaceAlpha = .2;
end 

if ~isfield(options,'EdgeCmap')
    options.EdgeCmap = hot(100);
end 

if ~isfield(options,'NodeCmap')
    options.NodeCmap = turbo(100);
end 

if ~isfield(options,'EdgeCmapLimits')
    options.EdgeCmapLimits = [min(A(:)) max(A(:))];
end 

if ~isfield(options,'NodeCmapLimits')
    options.NodeCmapLimits = [min(NodeValues) max(NodeValues)];
end 

if ~isfield(options,'PlotHemi')
    options.PlotHemi = 'Cortex';
end 

if ~isfield(options,'Surface')
    load('ICBMcortsurface.mat','cort_surface')
    options.Surface = cort_surface;
    clear cort_surface
end 

if ~isfield(options,'NodeMeasureLabel')
    options.NodeMeasureLabel = 'Nodal measure';
end 

if ~isfield(options,'EdgeMeasureLabel')
    options.EdgeMeasureLabel = 'Edge weight';
end

switch options.PlotHemi
    case 'Cortex'
        surface = options.Surface;
    case 'Left'
        VERTS = options.Surface.vertices;
        FACES = options.Surface.faces; 
        
        LVERTSIND = find(VERTS(:,1) < 0);
        
        LFACESIND = logical(sum(ismember(FACES,LVERTSIND),2) ~= 0);
        
        LVERTS = VERTS(LVERTSIND,:);
        LFACES = changem(FACES(LFACESIND,:),1:length(LVERTSIND),LVERTSIND);
        surface.vertices = LVERTS;
        surface.faces = LFACES;
        
        LCOORDS = coords(:,1) < 0;
        
        coords(~LCOORDS,:) = [];
        
        if size(NodeValues,2) == length(A)
            NodeValues(~LCOORDS) = [];
        end
        
        A(~LCOORDS,:) = [];
        A(:,~LCOORDS) = [];
        
    case 'Right'
        VERTS = options.Surface.vertices;
        FACES = options.Surface.faces; 
        
        RVERTSIND = find(VERTS(:,1) > 0);
        
        RFACESIND = logical(sum(ismember(FACES,RVERTSIND),2) ~= 0);
        
        RVERTS = VERTS(RVERTSIND,:);
        RFACES = changem(FACES(RFACESIND,:),1:length(RVERTSIND),RVERTSIND);
        surface.vertices = RVERTS;
        surface.faces = RFACES;
        
        RCOORDS = coords(:,1) > 0;
        
        coords(~RCOORDS,:) = [];
        
        if size(NodeValues,2) == length(A)
            NodeValues(~RCOORDS) = [];
        end
        
        A(~RCOORDS,:) = [];
        A(:,~RCOORDS) = [];   
        
end

p = patch(surface);
set(p,'EdgeColor','none','FaceColor',[.5 .5 .5],'FaceAlpha',options.SurfaceAlpha);
view([-90 0])
camlight;material dull

hold on
G = graph(A);
h = plot(G);
h.XData = coords(:,1);
h.YData = coords(:,2);
h.ZData = coords(:,3);

h.MarkerSize = options.NodeSize;
h.LineWidth = options.EdgeSize;
h.EdgeAlpha = options.EdgeAlpha;
h.NodeLabel = [];

if isempty(options.EdgeColor)
    if isfield(G.Edges,'Weight') || length(unique(A)) > 2
    h.EdgeColor = MapCmap(G.Edges.Weight,options.EdgeCmap,options.EdgeCmapLimits);
    else
       options.EdgeCmap = []; 
    end
    
else
    
    h.EdgeColor = options.EdgeColor;
    options.EdgeCmap = [];
    
end

if isempty(options.NodeColor)
    
    if length(NodeValues) == length(A)

        h.NodeColor = MapCmap(NodeValues,options.NodeCmap,options.NodeCmapLimits);
    
    else
        
        options.NodeCmap = [];
        
    end
    
else
    h.NodeColor = options.NodeColor;
    options.NodeCmap = [];
end

axis equal
axis off

cmap = [options.EdgeCmap; options.NodeCmap];

if ~isempty(cmap)

    colormap(cmap)

    NodeCmapN = size(options.NodeCmap,1);
    EdgeCmapN = size(options.EdgeCmap,1);
    cmapN = NodeCmapN + EdgeCmapN;

    EdgeMapProp = EdgeCmapN/cmapN;

    caxis([0 1])

    if NodeCmapN > 0
        
        NodeTickValues = ColorbarRange(NodeValues,options.NodeCmapLimits);
        
        c1 = colorbar('Limits',[EdgeMapProp 1]);

        c1.Ticks = rescale2range(NodeTickValues,EdgeMapProp,1,options.NodeCmapLimits(1),options.NodeCmapLimits(2));

        c1.TickLabels = NodeTickValues;

        c1.Label.String = options.NodeMeasureLabel;

    end

    if EdgeCmapN > 0 
        
        EdgeTickValues = ColorbarRange(G.Edges.Weight,options.EdgeCmapLimits);

        c2 = colorbar('Location','southoutside','Limits',[0 EdgeMapProp]);

        c2.Ticks = rescale2range(EdgeTickValues,0,EdgeMapProp,options.EdgeCmapLimits(1),options.EdgeCmapLimits(2));

        c2.TickLabels = EdgeTickValues;

        c2.Label.String = options.EdgeMeasureLabel;

    end

end

end

function Values = MapCmap(scalar,cmap,climits)
    
    if nargin < 3
       climits = [min(scalar) max(scalar)]; 
    end

    scalarclamped = scalar;
    scalarclamped(scalar < climits(1)) = climits(1);
    scalarclamped(scalar > climits(2)) = climits(2);
    Values = interp1(linspace(climits(1),climits(2),size(cmap,1)),cmap,scalarclamped);
          
end

function [TickValues] = ColorbarRange(data,Limits)

ax = axes;

imagesc(ax,data)

caxis(Limits)

c = colorbar;

TickValues = c.Ticks; 

delete(ax)

end

function B = rescale2range(A,l,u,inmin,inmax)

if nargin < 4

    inmin = min(A(:));

end

if nargin < 5

    inmax = max(A(:));

end

B = l + ((A-inmin)./(inmax-inmin)).*(u-l);

end