function PlotPropSurfSingle(Emp,MedLat,climits,cmap)

load('fsaverage_surface_data.mat','lh_inflated_verts','lh_faces','lh_rand200')

lhsurface.vertices = lh_inflated_verts;
lhsurface.faces = lh_faces;

%latpatch1 = axes('Position',[0.0100    0.5650    0.4315    0.4762]);
plotSurfaceROIBoundary(lhsurface,lh_rand200,Emp,'midpoint',cmap,1,2,climits);

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
