function Initalise_fetal_surface(BRAIN,hemi,POINTSPERJOB)

SURFACE_PARENT_DIR = '/projects/kg98/stuarto/Fetal_brain/CRL_Fetal_Brain_Atlas_2017';
OUTPUT_DIR = [SURFACE_PARENT_DIR,'/NewOptimFibreLength/',BRAIN];

mkdir(OUTPUT_DIR)

addpath(SURFACE_PARENT_DIR)
addpath([SURFACE_PARENT_DIR,'/gifti-1.6/'])

MSMFOLDER='MSM_outputs_smoothed_upsampled_hocr';
REGTYPE='sulc';

if strcmp(BRAIN,'Adult')
data = importorigasc([SURFACE_PARENT_DIR,'/Adult/surf/',hemi,'h.orig.asc']); 
nverts = data(1,1); data(1,:) = []; 
surface.vertices = single(data(1:nverts,1:3)); 
data(1:nverts,:) = []; 
surface.faces = data(:,1:3);     
  
else
fetal2adult = gifti([SURFACE_PARENT_DIR,'/',MSMFOLDER,'/',REGTYPE,'/MSM_',BRAIN,'_to_Adult_',hemi,'/target_deformed_',hemi,'_',BRAIN,'_to_Adult.surf.asc_anatresampled.surf.gii']); 
adult2fetal = gifti([SURFACE_PARENT_DIR,'/',MSMFOLDER,'/',REGTYPE,'/MSM_Adult_to_',BRAIN,'_',hemi,'/source_deformed_',hemi,'_Adult_to_',BRAIN,'.surf.asc_anatresampled.surf.gii']); 
nverts=size(fetal2adult.vertices,1);
v = zeros(nverts,3,2); 
v(:,:,1) = fetal2adult.vertices; 
v(:,:,2) = adult2fetal.vertices;

surface.vertices = mean(v,3); 
surface.faces = fetal2adult.faces;
    
end

[surf_redu,vertex_ind,surf_redu_points] = FibreDistanceSurface(surface,0.15,.1);

dlmwrite([OUTPUT_DIR,'/ReducedJobs.txt'],ceil(length(surf_redu.vertices) / POINTSPERJOB)); 

save([OUTPUT_DIR,'/Reduced.mat'],'surface','vertex_ind','surf_redu_points')
