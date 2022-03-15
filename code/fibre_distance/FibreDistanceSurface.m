function [surf_redu,vertex_ind,surf_redu_points] = FibreDistanceSurface(surface,reduce_factor,mm)

[surf_redu,vertex_ind] = reduce_surface(surface,reduce_factor);

surf_redu_points = PointsUnderSurface(surf_redu,faces,mm,1); 