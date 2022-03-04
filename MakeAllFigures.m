%% This script will produce all the figures in the paper

addpath(genpath('./'))

FIGURE_LOCATION = 'C:/Users/Stuart/Documents/GenerativeNetworkModel/Figures';

rand200_topo_mdls{1,1} = load('random200_FCV_TopoMdls_add3_Growth_0_output.mat');
rand200_topo_mdls{1,2} = load('random200_FCV_TopoMdls_add3_Growth_1_output.mat');
rand200_topo_mdls{2,1} = load('random200_FCV_TopoMdls_add2_Growth_0_output.mat');
rand200_topo_mdls{2,2} = load('random200_FCV_TopoMdls_add2_Growth_1_output.mat');
rand200_topo_mdls{3,1} = load('random200_FCV_TopoMdls_mult2_Growth_0_output.mat');
rand200_topo_mdls{3,2} = load('random200_FCV_TopoMdls_mult2_Growth_1_output.mat');

rand200_phys_mdls{2} = load('random200_FCV_PhysMdls_Growth_1_output.mat');
rand200_phys_mdls{1} = load('random200_FCV_PhysMdls_Growth_0_output.mat');

TOPO_MDL_LABELS = {'spatial','neighbours','matching','clu-avg','clu-min','clu-max','clu-diff','clu-prod','deg-avg','deg-min','deg-max','deg-diff','deg-prod'};
MDLTYPES = {'Spatial','Homophilly','Clustering','Degree'};

% First set of comparisons, each row gives the row and column of rand200_topo_mdls
% to pull data from
c1 = [1,2;2,1;1,1;1,1]; 

c2 = [1,1;3,1;2,1;3,1]; 

NCompareFigs = size(c1);

% There are two sets of every model i.e., growth and static, meaning there
% are 26 models. Comparing all combinations give 325 comparisons
NStatisticalComparisons = (26*25)/2;

FigureAnnot = {'','A','B','C'};

FigureNames = {'3','S1','S1','S1'};

GROUPNAMES = {{'Growth','Static'},{'Additive (no \gamma)','Multiplicative'},...
    {'Additive','Additive (no \gamma)'},{'Additive','Multiplicative'}};

MdlTypesInd = [1 2 2 3 3 3 3 3 4 4 4 4 4];

for i = 1:NCompareFigs

figure('Position',[233 364 1388 566])
ax1 = axes('Position',[0.0891    0.2915    0.9015    0.5566]);
[ordered_mdl_data{i},V,topo_mdl_order{i},issig,pval] = PlotMdlResults([rand200_topo_mdls{c1(i,1),c1(i,2)}.Fcv; rand200_topo_mdls{c2(i,1),c2(i,2)}.Fcv],TOPO_MDL_LABELS,'GrpNames',GROUPNAMES{i},...
'DataLabel','\it{F_{CV}}','SigLvl',NStatisticalComparisons,'MdlTypesInd',MdlTypesInd,'MdlTypesNames',MDLTYPES);

a = annotation('textbox',[0.0055    0.9468    0.0588    0.0885],'String',FigureAnnot{i},'EdgeColor','none','FontSize',48);

print([FIGURE_LOCATION,'/Figure',FigureNames{i},FigureAnnot{i},'.png'],'-dpng','-r300')

close all

end

%% Figure 4

LABELS = {'Spatial','Spatial+cCGE','cCGE','Spatial+uCGE','uCGE','MPC_{HIST}','Spatial+MPC_{HIST}','MPC_{T1/T2}','Spatial+MPC_{T1/T2}','Matching'};
PHYS_MDLTYPES = {'Spatial','Homophilly','cCGE','uCGE','MPC_{HIST}','MPC_{T1/T2}'};

phys_nStatisticalComparisons = (16*15)/2;

figure('Position',[233 364 1388 566])
ax1 = axes('Position',[0.0891    0.3888    0.9015    0.4593]);
GROWTH = rand200_phys_mdls{2}.Fcv;
STATIC = rand200_phys_mdls{1}.Fcv;
% The cCGE, uCGE, MPC(HIST) and MPC(T1T2) models actually ran in the growth
% form because I forgot they would. Anyway they are the exact same as the
% static models. Makes no sense to plot them, so I just remove them
for i = [3 5 6 8]
GROWTH{i} = [];
end
[phys_mdl_out,~,order_phys_mdls] = PlotMdlResults([GROWTH; STATIC],LABELS,'GrpNames',{'Growth','Static'},...
'DataLabel','\it{F_{CV}}','SigLvl',phys_nStatisticalComparisons,'MdlTypesInd',[1 3 3 4 4 5 5 6 6 2],'MdlTypesNames',PHYS_MDLTYPES);

ax1.XLabel.Position = [10.5000 -0.1865 1];

print([FIGURE_LOCATION,'/Figure4.png'],'-dpng','-r300')

close all

%% Figure 5
% This will take a very very long time to run.

% MdlColors = [152,78,163;...
% 55,126,184;...
% 228,26,28;...
% 77,175,74;...
% 255,127,0;...
% 66,224,245;...
% 166,86,40;...
% 77,190,238]./255;
% 
% load('random200_data4physmdl.mat')
% for i = 1:length(adjs)
% 
%     A = double(adjs{i}>0);
%     EmpData{1}(i,:) = sum(A); 
%     EmpData{3}(i,:) = betweenness_bin(A);
%     EmpData{2}(i,:) = clustering_coef_bu(A)';
%     EmpData{4}(i,:) = (sum(A_dist.*A))./EmpData{1}(i,:); 
% 
% end
% 
% rand200_phys_growth_topo = load('random200_FCV_PhysMdls_Growth_1_output_Topography.mat');
% 
% MDL = [4 5 10 1];
% 
% Model = rand200_phys_growth_topo.Topography{4};
% 
% CorrDATA = rand200_phys_mdls{2}.TopographyCorrs(MDL);
% 
% SurfaceData.vertices = lh_inflated_verts;
% SurfaceData.faces = lh_faces;
% SurfaceData.parc = lh_rand200;
% 
% load('Missing_rois.mat')
% SurfaceData.MissingROI = rand200_missing;
% 
% plotSpatialEmbeddings(EmpData,Model,CorrDATA,SurfaceData,LABELS(MDL),MdlColors([4 4 2 1],:),[FIGURE_LOCATION,'/Figure5']);
%
% close all

%% Figure S2

rand200_data = load('random200_data4topomdl.mat');
add3_optim_mdls = load('random200_TopoMdls_add3_Growth_0_output.mat', 'OptimMdls');
add2_optim_mdls = load('random200_TopoMdls_add2_Growth_0_output.mat', 'OptimMdls');
mult2_optim_mdls = load('random200_TopoMdls_mult2_Growth_0_output.mat', 'OptimMdls');

NETS{1} = mult2_optim_mdls.OptimMdls{3}.net;
NETS{2} = add2_optim_mdls.OptimMdls{3}.net;
NETS{3} = add3_optim_mdls.OptimMdls{3}.net;

for j = 1:3
for i = 1:length(NETS{j})
    a = zeros(100);
    a(NETS{j}{i}) = 1;
    a = a + a';
    NETS_FULL{j}{i} = a;
end
end

CompareEdgeLengtheCDF(rand200_data.A_dist,NETS_FULL{1},NETS_FULL{2},NETS_FULL{3},rand200_data.adjs)

exportgraphics(gcf,[FIGURE_LOCATION,'/FigureS2.png'],'resolution',300)

close all

%% Figure S3

VisualiseKSstats(rand200_topo_mdls{1,1}.KS,TOPO_MDL_LABELS,3,5)
print([FIGURE_LOCATION,'/FigureS3A.png'],'-dpng','-r300')

VisualiseKSstats(rand200_topo_mdls{1,2}.KS,TOPO_MDL_LABELS,3,5)
print([FIGURE_LOCATION,'/FigureS3B.png'],'-dpng','-r300')

close all

%% Figure S4

for i = 1:13
    
   alphavals_static_add3{i} = rand200_topo_mdls{1,1}.P{i}(:,3);
   alphavals_growth_add3{i} = rand200_topo_mdls{1,2}.P{i}(:,3);
   
   alphavals_static_add2{i} = rand200_topo_mdls{2,1}.P{i}(:,3);
   alphavals_growth_add2{i} = rand200_topo_mdls{2,2}.P{i}(:,3);
 
   etavals_static_add3{i} = rand200_topo_mdls{1,1}.P{i}(:,1);
   etavals_growth_add3{i} = rand200_topo_mdls{1,2}.P{i}(:,1); 
   
   gamvals_static_add3{i} = rand200_topo_mdls{1,1}.P{i}(:,2);
   gamvals_growth_add3{i} = rand200_topo_mdls{1,2}.P{i}(:,2);
 
   etavals_static_add2{i} = rand200_topo_mdls{2,1}.P{i}(:,1);
   etavals_growth_add2{i} = rand200_topo_mdls{2,2}.P{i}(:,1); 
   
   gamvals_static_mult2{i} = rand200_topo_mdls{3,1}.P{i}(:,2);
   gamvals_growth_mult2{i} = rand200_topo_mdls{3,2}.P{i}(:,2);
 
   etavals_static_mult2{i} = rand200_topo_mdls{3,1}.P{i}(:,1);
   etavals_growth_mult2{i} = rand200_topo_mdls{3,2}.P{i}(:,1);   
   
end

figure('Position',[233 364 1388 444])
ax1 = axes('Position',[0.0745    0.5286    0.92    0.4511]);
PlotMdlResults([etavals_growth_add3; etavals_static_add3],TOPO_MDL_LABELS,'GrpNames',{'Growth','Static'},...
'DataLabel','\eta','SigLvl',0,'MdlTypesInd',MdlTypesInd,'MdlTypesNames',MDLTYPES,'MdlOrder',topo_mdl_order{1});
a = annotation('textbox',[0.0055    0.9468    0.0588    0.0885],'String','A','EdgeColor','none','FontSize',48);

exportgraphics(gcf,[FIGURE_LOCATION,'/FigureS4A.png'],'resolution',300)

gamvals_static_add3{1} = [];
gamvals_growth_add3{1} = [];

figure('Position',[233 364 1388 444])
ax1 = axes('Position',[0.0745    0.5286    0.92    0.4511]);

PlotMdlResults([gamvals_growth_add3; gamvals_static_add3],TOPO_MDL_LABELS,'GrpNames',{'Growth','Static'},...
'DataLabel','\gamma','SigLvl',0,'MdlTypesInd',MdlTypesInd,'MdlTypesNames',MDLTYPES,'MdlOrder',topo_mdl_order{1});
a = annotation('textbox',[0.0055    0.9468    0.0588    0.0885],'String','B','EdgeColor','none','FontSize',48);

exportgraphics(gcf,[FIGURE_LOCATION,'/FigureS4B.png'],'resolution',300)

alphavals_static_add3{1} = [];
alphavals_growth_add3{1} = [];

figure('Position',[233 364 1388 444])
ax1 = axes('Position',[0.0745    0.5286    0.92    0.4511]);

PlotMdlResults([alphavals_growth_add3; alphavals_static_add3],TOPO_MDL_LABELS,'GrpNames',{'Growth','Static'},...
'DataLabel','\alpha','SigLvl',0,'MdlTypesInd',MdlTypesInd,'MdlTypesNames',MDLTYPES,'MdlOrder',topo_mdl_order{1});

a = annotation('textbox',[0.0055    0.9468    0.0588    0.0885],'String','C','EdgeColor','none','FontSize',48);

exportgraphics(gcf,[FIGURE_LOCATION,'/FigureS4C.png'],'resolution',300)

close all

%% Figure S5

scha200_phys_mdls{1} = load('Schaefer200_FCV_PhysMdls_Growth_0_output.mat');
scha200_phys_mdls{2} = load('Schaefer200_FCV_PhysMdls_Growth_1_output.mat');

LABELS = {'Spatial','Spatial+cCGE','cCGE','Spatial+uCGE','uCGE','MPC_{HIST}','Spatial+MPC_{HIST}','MPC_{T1/T2}','Spatial+MPC_{T1/T2}','Matching'};
PHYS_MDLTYPES = {'Spatial','Homophilly','cCGE','uCGE','MPC_{HIST}','MPC_{T1/T2}'};

phys_nStatisticalComparisons = (16*15)/2;

figure('Position',[233 364 1388 566])
ax1 = axes('Position',[0.0891    0.3888    0.9015    0.4593]);

GROWTH = scha200_phys_mdls{2}.Fcv;
STATIC = scha200_phys_mdls{1}.Fcv;
% The cCGE, uCGE, MPC(HIST) and MPC(T1T2) models actually ran in the growth
% form because I forgot they would. Anyway they are the exact same as the
% static models. Makes no sense to plot them, so I just remove them
for i = [3 5 6 8]
GROWTH{i} = [];
end
PlotMdlResults([GROWTH; STATIC],LABELS,'GrpNames',{'Growth','Static'},...
'DataLabel','\it{F_{CV}}','SigLvl',phys_nStatisticalComparisons,'MdlTypesInd',[1 3 3 4 4 5 5 6 6 2],'MdlTypesNames',PHYS_MDLTYPES);

ax1.XLabel.Position = [10.5000 -0.0812 1];

exportgraphics(gcf,[FIGURE_LOCATION,'/FigureS5.png'],'resolution',300)

close all

%% Figure S6

KS_PLOT_LABELS = LABELS;
for i = 1:10
     KS{i} = rand200_phys_mdls{1}.KS{i};
   KS_PLOT_LABELS{i} = ['Static ',LABELS{i}];
end

idx = 11;
for i = [1 2 4 7 9 10] 
         KS{idx} = rand200_phys_mdls{2}.KS{i};
   KS_PLOT_LABELS{idx} = ['Growth ',LABELS{i}];
   idx = idx + 1;
end
VisualiseKSstats(KS,KS_PLOT_LABELS,4,4)

exportgraphics(gcf,[FIGURE_LOCATION,'/FigureS6.png'],'resolution',300)

close all

%% Figure S7

rand200_phys_optim{1} = load('random200_PhysMdls_Growth_0_output.mat','DegCorr');
rand200_phys_optim{2} = load('random200_PhysMdls_Growth_1_output.mat','DegCorr');
figure('Position',[233 364 1388 536])
ax1 = axes('Position',[0.0891    0.3862    0.9015    0.4571]);
for i = 1:10
growth = [];
static = [];
for j = 1:100
growth = [growth; rand200_phys_optim{2}.DegCorr{i}{j}];
static = [static; rand200_phys_optim{1}.DegCorr{i}{j}];
end
GROWTH_DEGCORR{i} = growth;
STATIC_DEGCORR{i} = static;
end
for i = [3 5 6 8]
GROWTH_DEGCORR{i} = [];
end

PlotMdlResults([GROWTH_DEGCORR; STATIC_DEGCORR],LABELS,'GrpNames',{'Growth','Static'},...
'DataLabel',{'Spearman correlation','with empirical degree'},'SigLvl',0,'MdlTypesInd',[1 3 3 4 4 5 5 6 6 2],'MdlTypesNames',PHYS_MDLTYPES,'MdlOrder',order_phys_mdls);

ax1.XLabel.Position = [10.5000   -1.7843    1.0000];

exportgraphics(gcf,[FIGURE_LOCATION,'/FigureS7.png'],'resolution',300)

close all