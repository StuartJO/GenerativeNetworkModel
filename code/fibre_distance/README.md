# Making fibre distances

This code is calculates fibre distances as proposed in Modeling spatial, developmental, physiological, and topological constraints on human brain connectivity

To start with this is how I ran these scripts (largely for my own reference)

```
% In MATLAB
TIMEPOINTS = {'21','22','23','24','25','26','27','28','29','30','31','32','33','34','35','36','37','38','Adult'};
for i = 1:length(TIMEPOINTS)
    Initalise_fetal_surface(TIMEPOINTS{i},'l',100)
end
```

Then in the terminal run
```
POINTSPERJOB=100
for BRAIN in {21..38} Adult; do
DATADIR="/projects/kg98/stuarto/Fetal_brain/CRL_Fetal_Brain_Atlas_2017/NewOptimFibreLength/${BRAIN}"
CODEDIR="/fs02/hf49/Stuart/GrowthModel_newParc/GenerativeNetworkModel"
${CODEDIR}/code/fibre_distance/RunFindDirectConnections.sh ${DATADIR} ${CODEDIR} ${POINTSPERJOB}
done
```
Ensure the above all ran successfully then run

```
NODESPERJOB=100
for BRAIN in {21..38} Adult; do
DATADIR="/projects/kg98/stuarto/Fetal_brain/CRL_Fetal_Brain_Atlas_2017/NewOptimFibreLength/${BRAIN}"
CODEDIR="/fs02/hf49/Stuart/GrowthModel_newParc/GenerativeNetworkModel"
${CODEDIR}/code/fibre_distance/RunFindShortestPath.sh ${DATADIR} ${CODEDIR} ${NODESPERJOB}
done
```

And finally run in MATLAB
```
TIMEPOINTS = {'21','22','23','24','25','26','27','28','29','30','31','32','33','34','35','36','37','38','Adult'};
for i = 1:length(TIMEPOINTS)
    DATADIR='[/projects/kg98/stuarto/Fetal_brain/CRL_Fetal_Brain_Atlas_2017/NewOptimFibreLength/',TIMEPOINTS{i}];
    dists{i} = GetParcellationDistance(DATADIR,lh_rand200)
end
```

Now if you want to perform this on your own data, you will need to adapt the above. Luckily for you, I have made a headstart:

```
% Need to face the faces and vertices for your surface
surface.faces = faces;
surface.vertices = vertices;

% Define your data directory (where all outputs for fibre_distance scripts will output
DATADIR=''

% Reduce the surface to retain only 15% of vertices, and find points .1mm under the surface
[surf_redu,vertex_ind,surf_redu_points] = FibreDistanceSurface(surface,0.15,.1);

% To speed things up (because by god this script desperately needs it), we 
% divide up calculations so they are performed across multiple slurm jobs.
% The calculations are performed over the vertices on the reduced surface 
% (surf_redu). I found performing calculations over 100 vertices per job 
% offered "good" enough performance

POINTSPERJOB=100;
TOTALJOBS=ceil(length(surf_redu.vertices) / POINTSPERJOB)

dlmwrite([DATADIR,'/ReducedJobs.txt'],TOTALJOBS); 

save([OUTPUT_DIR,'/Reduced.mat'],'surface','vertex_ind','surf_redu_points')
```
Then you will need to run the following scripts in this order

```
POINTSPERJOB=100
DATADIR="YOUR_DATA_LOCATION_HERE"
CODEDIR="CODE_LOCATION_HERE"
${CODEDIR}/code/fibre_distance/RunFindDirectConnections.sh ${DATADIR} ${CODEDIR} ${POINTSPERJOB}
```
Ensure the above all ran successfully then run

```
NODESPERJOB=100
DATADIR="YOUR_DATA_LOCATION_HERE"
CODEDIR="CODE_LOCATION_HERE"
${CODEDIR}/code/fibre_distance/RunFindShortestPath.sh ${DATADIR} ${CODEDIR} ${NODESPERJOB}
```

And finally run in MATLAB
```
%parcdata needs to be the ROI id for each vertex in the original surface    
parcdata=your_parc_data;
DATADIR='YOUR_DATA_LOCATION_HERE';
dists = GetParcellationDistance(DATADIR,parcdata)
```