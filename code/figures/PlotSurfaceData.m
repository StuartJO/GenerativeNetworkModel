function PlotSurfaceData(Data,MedLat,surface,climits,cmap)

plotSurfaceROIBoundary(surface,surface.parc,Data,'midpoint',cmap,1,2,climits);

if MedLat == 1
camlight(80,-10);
camlight(-80,-10);
view([-90 0])
axis off
axis tight
axis equal

else
    
camlight(80,-10);
camlight(-80,-10);
view([90 0])
axis off
axis tight
axis equal

end
