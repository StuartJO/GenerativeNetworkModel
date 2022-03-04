# GenerativeNetworkModel

This code is for Modeling spatial, developmental, physiological, and topological constraints on human brain connectivity

Any questions please email stuart.oldham@mcri.edu.au

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

You may want to use this code for your own data. Use the runGenPhysMdl.m as a template, as I used this as a wrapped for all my topological model. Essentially all you need is a) a number of connectomes for a number of subjects arranged in a cell (adjs), a final distance matrix to evaluate the model with (A_dist, code currently assumes only a single A_dist applies to all networks but you could have individualised ones)), and finally a distance matrix to run the model with (if not using a growth model, this will be the same as A_dist). It should be very straightforward to adjust the code for your own purpose! 


