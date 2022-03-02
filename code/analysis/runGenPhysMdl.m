function runGenPhysMdl(SUB,GROWTH,PARC,SAVEDIR,MDLS)

if nargin < 4
   SAVEDIR = './'; 
end

if nargin < 5
    MDLS = 1:10;
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
    D = mdldata.dists;
else
    D = A_dist;
end   
    
MDLIND = 1;

for MDL = MDLS

if MDL == 1
display(num2str(MDL))
Input.PD1Func = 'exponential';
Input.PD2Func = 'powerlaw';
Input.ModelNum=1;
Input.AddMult = 'Mult';
PD = [];
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

PD = mdldata.cCGE;

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

D = mdldata.cCGE;
PD = [];

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

PD = mdldata.uCGE;

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

D = mdldata.uCGE;
PD = [];
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
    

D = mdldata.hist_mpc;
PD = [];
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
    
PD = mdldata.hist_mpc;

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
    

D = mdldata.t1t2_mpc;
PD = [];
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
    
PD = mdldata.t1t2_mpc;
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

Output = GenMdl(A,A_dist,D,PD,Input);

if length(MDLS) > 1
    Outputs.maxKS{MDLIND} = Output.maxKS;
    Outputs.KS{MDLIND} = Output.KS;
    Outputs.P{MDLIND} = Output.P;
    Outputs.b{MDLIND} = Output.b;
    Outputs.maxKS_optim{MDLIND} = Output.maxKS_optim;
    Outputs.KS_optim{MDLIND} = Output.KS_optim;
    Outputs.b_optim{MDLIND} = Output.b_optim;

% Save the input configurations to output. Helps to keep track of what was
% done
Outputs.Input{MDLIND} = Output.Input;
    MDLIND = MDLIND + 1;
end

end

if length(MDLS) > 1
    save([SAVEDIR,parcname,'_PhysMdls_Sub_',num2str(SUB),'_add3_Growth_',num2str(GROWTH),'.mat'],'-struct','Outputs','-v7.3')
else
    save([SAVEDIR,parcname,'_PhysMdls_Mdl_',num2str(MDL),'_Sub_',num2str(SUB),'_add3_Growth_',num2str(GROWTH),'.mat'],'-struct','Output','-v7.3')
end
