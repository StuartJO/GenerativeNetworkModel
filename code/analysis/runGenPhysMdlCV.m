function runGenPhysMdlCV(PARC,GROWTH,MdlNum,ITER)

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

display(['Performing cross-validation for ',parcname,', phys growth ',num2str(GROWTH),' model number ',num2str(MdlNum),', iteration ',num2str(ITER)]);

adjs = mdldata.adjs;
A_dist = mdldata.A_dist;
    
if GROWTH
    D = mdldata.dists;
else
    D = A_dist;
end   

% Set up some default options
% Input.Growth = GROWTH;
% Input.PD1Func = 'exponential';
% Input.IndvDist = 0;
% Input.PD2Func = 'powerlaw';
% Input.epsilon = 0;
% Input.seed = [];
% Input.normsum = 0;
% Input.TopoFunc = 'powerlaw';

if MdlNum == 1

% Input.PD1Func = 'exponential';
% Input.PD2Func = 'powerlaw';
% Input.ModelNum=1;
% Input.AddMult = 'Mult';
PD = [];

elseif MdlNum == 2
     
% Input.PD1Func = 'exponential';
% Input.PD2Func = 'powerlaw';
% Input.ModelNum=1;
% Input.AddMult = 'Add';

PD = mdldata.cCGE;

elseif MdlNum == 3

% Input.ModelNum=1;
% Input.AddMult = 'Mult';
% Input.PD1Func = 'powerlaw';
% Input.PD2Func = 'powerlaw';

D = mdldata.cCGE;
PD = [];

elseif MdlNum == 4

% Input.ModelNum=1;
% Input.AddMult = 'Add';
% Input.PD1Func = 'exponential';
% Input.PD2Func = 'powerlaw';

PD = mdldata.uCGE;

elseif MdlNum == 5

% Input.ModelNum=1;
% Input.AddMult = 'Mult';
% Input.PD1Func = 'powerlaw';
% Input.PD2Func = 'powerlaw';

D = mdldata.uCGE;
PD = [];

elseif MdlNum == 6
    
% Input.PD1Func = 'powerlaw';
% Input.PD2Func = 'powerlaw';
% Input.AddMult = 'Mult';
% Input.ModelNum=1;

D = mdldata.hist_mpc;
PD = [];


elseif MdlNum == 7

% Input.PD1Func = 'exponential';
% Input.PD2Func = 'powerlaw';
% Input.ModelNum=1;
% Input.AddMult = 'Add';
    
PD = mdldata.hist_mpc;

elseif MdlNum == 8
    
% Input.PD1Func = 'powerlaw';
% Input.PD2Func = 'powerlaw';
% Input.ModelNum=1;
% Input.AddMult = 'Mult';

D = mdldata.t1t2_mpc;
PD = [];


elseif MdlNum == 9

% Input.PD1Func = 'exponential';
% Input.PD2Func = 'powerlaw';
% Input.ModelNum=1;
% Input.AddMult = 'Add';
    
PD = mdldata.t1t2_mpc;


elseif MdlNum == 10

% Input.AddMult = 'Add';
% Input.ModelNum=3;
PD = [];

end

Mdlouts = load([parcname,'_PhysMdls_Growth_',num2str(GROWTH),'_output.mat']);

Input = Mdlouts.Inputs{MdlNum};

P = Mdlouts.OptimMdl{MdlNum}.min_maxKS.P;
Fcv = CrossValidateModel(adjs,A_dist,D,PD,P,1,Input);

Fcv.P = P;

output_location = '/fs02/hf49/Stuart/GrowthModel_newParc/Paper_schaefer_mdls/Crossvalidated/';

save([output_location,'CrossValidate_',parcname,'_PhysMdls_mdl_',num2str(MdlNum),'_Growth_',num2str(GROWTH),'_iter_',num2str(ITER),'.mat'],'-struct','Fcv','-v7.3')
