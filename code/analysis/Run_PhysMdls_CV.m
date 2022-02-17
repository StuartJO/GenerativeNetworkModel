function Run_PhysMdls_CV(PARC,GROWTH,MdlNum,ITER)

Input.Growth = GROWTH;
Input.DistFunc = 'exponential';
Input.IndvDist = 0;
Input.GeneFunc = 'powerlaw';

if PARC == 1
mdldata = load('Schaefer200_data4physmdl.mat');
parcname = 'Schaefer200';
elseif PARC == 2
mdldata = load('Schaefer400_data4physmdl.mat');
parcname = 'Schaefer400';
elseif PARC == 3
mdldata = load('random200_data4physmdl.mat');
parcname = 'random200';
elseif PARC == 4
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
if MdlNum == 1

Input.DistFunc = 'exponential';
Input.GeneFunc = 'powerlaw';
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

elseif MdlNum == 2
 
% Additive 3 Spatial + corrected gene         
Input.DistFunc = 'exponential';
Input.GeneFunc = 'powerlaw';
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
elseif MdlNum == 3

Input.ModelNum=1;
Input.AddMult = 'Mult';
Input.DistFunc = 'powerlaw';
Input.GeneFunc = 'powerlaw';

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
elseif MdlNum == 4


Input.ModelNum=1;
Input.AddMult = 'Add';
Input.DistFunc = 'exponential';
Input.GeneFunc = 'powerlaw';

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
elseif MdlNum == 5


Input.ModelNum=1;
Input.AddMult = 'Mult';
Input.DistFunc = 'powerlaw';
Input.GeneFunc = 'powerlaw';

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

elseif MdlNum == 6
    
Input.DistFunc = 'powerlaw';
Input.GeneFunc = 'powerlaw';
Input.normsum = 0;
Input.ParamRange(1,:) = [0 50];
Input.ModelNum=1;
    

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

elseif MdlNum == 7

Input.DistFunc = 'exponential';
Input.GeneFunc = 'powerlaw';
Input.normsum = 0;
Input.ParamRange(1,:) = [-2 0];
Input.ModelNum=1;
    
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


elseif MdlNum == 8
    
Input.DistFunc = 'powerlaw';
Input.GeneFunc = 'powerlaw';
Input.normsum = 0;
Input.ParamRange(1,:) = [0 50];
Input.ModelNum=1;
    

D = mdldata.t1t2_mpc;
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

elseif MdlNum == 9

Input.DistFunc = 'exponential';
Input.GeneFunc = 'powerlaw';
Input.normsum = 0;
Input.ParamRange(1,:) = [-2 0];
Input.ModelNum=1;
    
PD = mdldata.t1t2_mpc;
Input.AddMult = 'Add';
%fileoutname = 'add3';
% gamma
% gamma
Input.ParamRange(2,:) = [0 0];
% lambda (PD2 param)
Input.ParamRange(5,:) = [0 50];
% alpha
Input.ParamRange(3,:) = [0 0];
% alpha2 (alpha for PD2)
Input.ParamRange(4,:) = [0 8];  

elseif MdlNum == 10

Input.AddMult = 'Add';
Input.ModelNum=3;
PD = [];

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

Mdlouts = load([parcname,'_PhysMdls_Growth_',num2str(GROWTH),'_output.mat']);

P = Mdlouts.MinMdlFitParams{MdlNum};
Fcv = CrossValidateModel(adjs,A_dist,D,PD,P,1,Input);

Fcv.P = P;

output_location = '/fs02/hf49/Stuart/GrowthModel_newParc/Paper_schaefer_mdls/Crossvalidated/';

save([output_location,'CrossValidate_',parcname,'_PhysMdls_mdl_',num2str(MdlNum),'_Growth_',num2str(GROWTH),'_iter_',num2str(ITER),'.mat'],'-struct','Fcv','-v7.3')
