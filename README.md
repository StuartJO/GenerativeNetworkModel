# GenerativeNetworkModel

This code is for Modeling spatial, developmental, physiological, and topological constraints on human brain connectivity

This requires the Brain Connectivity Toolbox (https://sites.google.com/site/bctnet/) to run some functions (place it into this directory).

Any questions please email stuart.oldham@mcri.edu.au

Data to regenerate all figures and data reported in the paper can be found here:
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6341625.svg)](https://doi.org/10.5281/zenodo.6341625)

Note that for the posted data, the raw optimisation and cross validation outputs are not available because it is too much data to host anywhere online. The optimisation data runs to about 120GB and the cross validation to 400GB, which is a lot! If you you really, *really*, ***really*** want to access it, chat to me and we will see what we can figure out.

To generative all figures from the paper run (in MATLAB)

```
MakeAllFigures.m
```

To rerun all analysis for the paper, run (from a terminal)

```
mkdir ./data/Optimisation
for TYPE in 1 2 3; do
    for GROWTH in 0 1; do
        sbatch ./code/RunTopoMdls.sh $GROWTH $TYPE
    done
done

for PARC in 1 2; do
    for GROWTH in 0 1; do
        sbatch ./code/RunPhysMdls.sh $GROWTH $PARC
    done
done
```

Wait 4-5 weeks, then run in MATLAB


```
Compile_TopoMdls
Compile_PhysMdls
```

Then after that run (from a terminal)

```
mkdir ./data/Crossvalidated
for TYPE in 1 2 3; do
    for GROWTH in 0 1; do
        sbatch ./code/RunTopoMdlsCV.sh $GROWTH $TYPE
    done
done

for PARC in 1 2; do
    for GROWTH in 0 1; do
        sbatch ./RunPhysMdlsCV.sh $GROWTH $PARC
    done
done
```

Wait another 1-2 weeks and then finally run (from MATLAB)

```
Compile_TopoMdls_FCV
Compile_PhysMdls_FCV
```

et voilà! You have rerun all the analysis! Yayyyy

You may want to use this code for your own data. I have provided some functions which allow you to do so. For topological models, you can do the following:

```
% A = cell of adjacency networks to model
% A_dist = a final distance matrix to use to calculate maxKS
% D = a distance matrix to use in the model itself, or a cell array where 
% each cell contains a distance matrix. Probably for your data this is the same as A_dist
addpath(genpath(('./')))

% If you want to have an individualised seed network for each subject, the code isn't set up to do that (sorry not sorry)

Nsubs = length(A);
mkdir ./data/CustomOptimisation/
for sub = 1:Nsubs
MdlOutput = runCustomGenTopoMdl(A{sub},A_dist,D,'add3',1:13);
save(['./data/CustomOptimisation/Sub_',num2str(i),'_add3_topomdls.mat'],'-v7.3','struct','MdlOutput')
end

%% Run cross validation

fileformat = ['./data/CustomOptimisation/Sub_#_add3_topomdls.mat'];
        
CompiledOutput = CompileGenMdlOutputs(fileformat,1:Nsubs,1);

mkdir ./data/CustomCrossValidation/

for mdl = 1:13
    for iter = 1:20
    P = CompiledOutput.OptimMdl{mdl}.min_maxKS.P;
    Input = CompiledOutput.Inputs{mdl};
    CV_output = CrossValidateModel(A,A_dist,D,[],P,1,Input);
    save(['./data/CustomCrossValidation/Mdl_',num2str(mdl),'_add3_topomdls_iter_',num2str(iter),'.mat'],'-v7.3','struct','MdlOutput')
    end
end
fileformat = ['./data/CustomCrossValidation/Mdl_#_add3_topomdls_iter_$.mat'];
CV = CompileCVOutputs(fileformat,1:13,1,Nsubs,0);
```
I want to stress the above is purely an illustrative example. Running the above exactly as shown in MATLAB will take a very very very long time. You can easily adapt the above to be seperate functions which you can then run on a cluster (ADVISED).

If you have same other measure of node/region similarity/distance and want to see how that does (much like we had distance and gene expression for example) you can try the following:

```
% A = cell of adjacency networks to model
% A_dist = a final distance matrix to use to calculate maxKS
% PD1 = a distance matrix to use in the model itself, or a cell array where 
% each cell contains a distance matrix. Probably for your data this is the same as A_dist
% PD2 = a distance/similarity matrix to use in the model itself, or a cell array where 
% each cell contains a distance/similarity matrix.


% These are the same parameters I use to run the spatial+uCGE model
Input.ModelNum=1;
% Model formulation
Input.AddMult = 'Add';

Input.PD1Func = 'exponential';
Input.PD2Func = 'powerlaw';

% eta (PD1 param)
Input.ParamRange(1,:) = [-2 0];

% gamma
Input.ParamRange(2,:) = [NaN NaN];

% lambda (PD2 param)
Input.ParamRange(5,:) = [-50 250];

% alpha (for topology)
Input.ParamRange(3,:) = [NaN NaN];

% alpha2 (alpha for PD2)
Input.ParamRange(4,:) = [0 8];

addpath(genpath(('./')))

Nsubs = length(A);
mkdir ./data/CustomOptimisation/
for sub = 1:Nsubs
MdlOutput = runCustomGenPhysMdl(A{sub},A_dist,PD1,PD2,Input);
save(['./data/CustomOptimisation/Sub_',num2str(i),'_physmdl.mat'],'-v7.3','struct','MdlOutput')
end

%% Run cross validation

fileformat = ['./data/CustomOptimisation/Sub_#_physmdl.mat'];
        
CompiledOutput = CompileGenMdlOutputs(fileformat,1:Nsubs,1);

mkdir ./data/CustomCrossValidation/

for mdl = 1
    for iter = 1:20
    P = CompiledOutput.OptimMdl{mdl}.min_maxKS.P;
    Input = CompiledOutput.Inputs{mdl};
    CV_output = CrossValidateModel(A,A_dist,PD1,PD2,P,1,Input);
    save(['./data/CustomCrossValidation/Mdl_',num2str(mdl),'_physmdls_iter_',num2str(iter),'.mat'],'-v7.3','struct','MdlOutput')
    end
end
fileformat = ['./data/CustomCrossValidation/Mdl_#_physmdls_iter_$.mat'];
CV = CompileCVOutputs(fileformat,1,1,Nsubs,0);

```
As before, this is just illustrative! Advise writing the above into functions you can submit on a cluster
