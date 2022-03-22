function runGenPhysMdl(SUB,GROWTH,PARC,SAVELOC,MDLS,save_suffix)

% This is a wrapper for all the physiological models in the paper.
%
% Inputs:
%
% SUB = subject to run the models on
%
% GROWTH = 1 to run a growth model, 0 for a static model
%
% PARC = value between 1-3 to indicate which parcellation/data to use,1 =
%   random200, 2 = Schaefer200, 3 = random200 with 18 timepoints (remove
%   adult timepoints)
%
% SAVELOC = location to save output (optional, will save in current 
%   directory by default)
%
% MDLS = index of models to run (optional, defaults to 1:10)
%
% save_suffix = a string to add to the end of the output filename

if nargin < 4
   SAVELOC = '.'; 
end

if nargin < 5
    MDLS = 1:10;
end

if nargin < 6
   save_suffix = ''; 
end

Input.Growth = GROWTH;
Input.PD1Func = 'exponential';
Input.PD2Func = 'powerlaw';
Input.useParfor = 1;

Input.TopoFunc = 'powerlaw';
Input.seed = [];

Input.normsum = 0;
Input.ndraw = 2000;
Input.pow = 2;
Input.nlvl = 5;

PHYS_MDL_NAME = {'Spatial','Spatial+cCGE','cCGE','Spatial+uCGE','uCGE','MPC_{HIST}','Spatial+MPC_{HIST}','MPC_{T1/T2}','Spatial+MPC_{T1/T2}','Matching'};

if PARC == 1
mdldata = load('random200_data4physmdl.mat');
parcname = 'random200';
elseif PARC == 2
mdldata = load('Schaefer200_data4physmdl.mat');
parcname = 'Schaefer200';
elseif PARC == 3
mdldata = load('random200_data4physmdl_18timepoints.mat');
parcname = 'random200_18tps';    
end

A = double(mdldata.adjs{SUB}>0);

A_dist = mdldata.A_dist;
    
if GROWTH
    PD1 = mdldata.dists;
else
    PD1 = A_dist;
end   
    
MDLIND = 1;

for MDL = MDLS

if MDL == 1
display(num2str(MDL))
Input.PD1Func = 'exponential';
Input.PD2Func = 'powerlaw';
Input.ModelNum=1;
% Note that for a single parameter model, additive and multiplicative are
% the same. The additive formulation does an extra bit of computation in
% normalising PD1, thus taking slightly more time, but this normalisation
% has no effect when no other term is resolved. Due to it taking slightly
% longer (and because I had many repeats to do) I just use the
% multiplicative code for single parameter models
Input.AddMult = 'Mult';
PD2 = [];
% eta (PD1 param)
Input.ParamRange(1,:) = [-2 0];
% gamma
Input.ParamRange(2,:) = [0 0];
% lambda (PD2 param)
Input.ParamRange(5,:) = [0 0];
% alpha
Input.ParamRange(3,:) = [1 1];
% alpha2 (alpha for PD2)
Input.ParamRange(4,:) = [1 1];

elseif MDL == 2
 
% Additive 3 Spatial + corrected gene         
Input.PD1Func = 'exponential';
Input.PD2Func = 'powerlaw';
display(num2str(MDL))
Input.ModelNum=1;
Input.AddMult = 'Add';

PD2 = mdldata.cCGE;

% eta (PD1 param)
Input.ParamRange(1,:) = [-2 0];
% gamma
Input.ParamRange(2,:) = [1 1];
% lambda (PD2 param)
Input.ParamRange(5,:) = [-50 250];
% alpha
Input.ParamRange(3,:) = [0 0];
% alpha2 (alpha for PD2)
Input.ParamRange(4,:) = [0 8];

 % Corrected gene
elseif MDL == 3
display(num2str(MDL))

Input.ModelNum=1;
Input.AddMult = 'Mult';
Input.PD1Func = 'powerlaw';
Input.PD2Func = 'powerlaw';

PD1 = mdldata.cCGE;
PD2 = [];

% eta (PD1 param)
Input.ParamRange(1,:) = [-50 100];
% gamma
Input.ParamRange(2,:) = [0 0];
% lambda (PD2 param)
Input.ParamRange(5,:) = [0 0];
% alpha
Input.ParamRange(3,:) = [0 0];
% alpha2 (alpha for PD2)
Input.ParamRange(4,:) = [0 0];

% Additive 3 Spatial + non corrected gene 
elseif MDL == 4
display(num2str(MDL))

Input.ModelNum=1;
Input.AddMult = 'Add';
Input.PD1Func = 'exponential';
Input.PD2Func = 'powerlaw';

PD2 = mdldata.uCGE;

% eta (PD1 param)
Input.ParamRange(1,:) = [-2 0];
% gamma
Input.ParamRange(2,:) = [1 1];
% lambda (PD2 param)
Input.ParamRange(5,:) = [-50 250];
% alpha
Input.ParamRange(3,:) = [0 0];
% alpha2 (alpha for PD2)
Input.ParamRange(4,:) = [0 8];
% Non Corrected gene
elseif MDL == 5
display(num2str(MDL))

Input.ModelNum=1;
Input.AddMult = 'Mult';
Input.PD1Func = 'powerlaw';
Input.PD2Func = 'powerlaw';

PD1 = mdldata.uCGE;
PD2 = [];
% eta (PD1 param)
Input.ParamRange(1,:) = [-50 100];
% gamma
Input.ParamRange(2,:) = [0 0];
% lambda (PD2 param)
Input.ParamRange(5,:) = [0 0];
% alpha
Input.ParamRange(3,:) = [0 0];
% alpha2 (alpha for PD2)
Input.ParamRange(4,:) = [0 0];

elseif MDL == 6
    
Input.PD1Func = 'powerlaw';
Input.PD2Func = 'powerlaw';
Input.normsum = 0;
Input.ParamRange(1,:) = [0 50];
    

PD1 = mdldata.hist_mpc;
PD2 = [];
      Input.AddMult = 'Mult';
%    fileoutname = 'mult2';
% gamma
Input.ParamRange(2,:) = [0 0];
% lambda (PD2 param)
Input.ParamRange(5,:) = [0 0];
% alpha
Input.ParamRange(3,:) = [1 1];
% alpha2 (alpha for PD2)
Input.ParamRange(4,:) = [1 1];   

elseif MDL == 7

Input.PD1Func = 'exponential';
Input.PD2Func = 'powerlaw';
Input.normsum = 0;
Input.ParamRange(1,:) = [-2 0];
    
PD2 = mdldata.hist_mpc;

Input.AddMult = 'Add';

% gamma
Input.ParamRange(2,:) = [0 0];
% lambda (PD2 param)
Input.ParamRange(5,:) = [0 50];
% alpha
Input.ParamRange(3,:) = [0 0];
% alpha2 (alpha for PD2)
Input.ParamRange(4,:) = [0 8];   


elseif MDL == 8
    
Input.PD1Func = 'powerlaw';
Input.PD2Func = 'powerlaw';
Input.normsum = 0;
Input.ParamRange(1,:) = [0 50];
    

PD1 = mdldata.t1t2_mpc;
PD2 = [];
Input.AddMult = 'Mult';
% gamma
Input.ParamRange(2,:) = [0 0];
% lambda (PD2 param)
Input.ParamRange(5,:) = [0 0];
% alpha
Input.ParamRange(3,:) = [1 1];
% alpha2 (alpha for PD2)
Input.ParamRange(4,:) = [1 1];   

elseif MDL == 9

Input.ModelNum=1;
Input.PD1Func = 'exponential';
Input.PD2Func = 'powerlaw';
Input.normsum = 0;
Input.ParamRange(1,:) = [-2 0];
    
PD2 = mdldata.t1t2_mpc;
Input.AddMult = 'Add';
% gamma
Input.ParamRange(2,:) = [0 0];
% lambda (PD2 param)
Input.ParamRange(5,:) = [0 50];
% alpha
Input.ParamRange(3,:) = [0 0];
% alpha2 (alpha for PD2)
Input.ParamRange(4,:) = [0 8];  

elseif MDL == 10

Input.ModelNum=3;
Input.AddMult = 'Add';
Input.PD1Func = 'exponential';
Input.TopoFunc = 'powerlaw';
PD2 = [];
% eta (PD1 param)
Input.ParamRange(1,:) = [-2 0];
% gamma
Input.ParamRange(2,:) = [-8 8];
% lambda (PD2 param)
Input.ParamRange(5,:) = [0 0];
% alpha
Input.ParamRange(3,:) = [0 8];
% alpha2 (alpha for PD2)
Input.ParamRange(4,:) = [0 0];   

end

Output = GenMdl(A,A_dist,PD1,PD2,Input);
Output.Input.PhysMdlName = PHYS_MDL_NAME{MDL};

Outputs.maxKS{MDLIND} = Output.maxKS;
Outputs.DegCorr{MDLIND} = Output.DegCorr;
Outputs.KS{MDLIND} = Output.KS;
Outputs.P{MDLIND} = Output.P;
Outputs.b{MDLIND} = Output.b;

Outputs.optim_maxKS{MDLIND} = Output.optim_maxKS;
Outputs.optim_KS{MDLIND} = Output.optim_KS;
Outputs.optim_b{MDLIND} = Output.optim_b;
Outputs.optim_DegCorr{MDLIND} = Output.optim_DegCorr;

Outputs.bestDegCorr_maxKS{MDLIND} = Output.bestDegCorr_maxKS;
Outputs.bestDegCorr_KS{MDLIND} = Output.bestDegCorr_KS;
Outputs.bestDegCorr_b{MDLIND} = Output.bestDegCorr_b;
Outputs.bestDegCorr_DegCorr{MDLIND} = Output.bestDegCorr_DegCorr;

% Save the input configurations to output. Helps to keep track of what was
% done
Outputs.Input{MDLIND} = Output.Input;
    MDLIND = MDLIND + 1;

end

save([SAVELOC,'/',parcname,'_PhysMdls_Sub_',num2str(SUB),'_add3_Growth_',num2str(GROWTH),save_suffix,'.mat'],'-struct','Outputs','-v7.3')

