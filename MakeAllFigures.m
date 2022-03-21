%% This script will produce all the figures in the paper

addpath(genpath('./'))

FIGURE_LOCATION = 'C:/Users/Stuart/Documents/GenerativeNetworkModel/Figures';

rand200_topo_mdls{1,1} = load('random200_FCV_TopoMdls_add3_Growth_0_output.mat');
rand200_topo_mdls{1,2} = load('random200_FCV_TopoMdls_add3_Growth_1_output.mat');
rand200_topo_mdls{2,1} = load('random200_FCV_TopoMdls_add2_Growth_0_output.mat');
rand200_topo_mdls{2,2} = load('random200_FCV_TopoMdls_add2_Growth_1_output.mat');
rand200_topo_mdls{3,1} = load('random200_FCV_TopoMdls_mult2_Growth_0_output.mat');
%rand200_topo_mdls{3,2} = load('random200_FCV_TopoMdls_mult2_Growth_1_output.mat');

rand200_phys_mdls{2} = load('random200_FCV_PhysMdls_Growth_1_output.mat');
rand200_phys_mdls{1} = load('random200_FCV_PhysMdls_Growth_0_output.mat');

TOPO_MDL_LABELS = {'spatial','neighbours','matching','clu-avg','clu-min','clu-max','clu-diff','clu-prod','deg-avg','deg-min','deg-max','deg-diff','deg-prod'};
MDLTYPES = {'Spatial','Homophilly','Clustering','Degree'};

% First set of comparisons, each row gives the row and column of rand200_topo_mdls
% to pull data from
c1 = [1,2;2,1;1,1;1,1;1,2]; 

c2 = [1,1;3,1;2,1;3,1;2,2]; 

NCompareFigs = size(c1);

% There are two sets of every model i.e., growth and static, meaning there
% are 26 models. Comparing all combinations give 325 comparisons
NStatisticalComparisons = (26*25)/2;

FigureAnnot = {'','A','B','C',''};

FigureNames = {'3','S1','S1','S1','S3'};

GROUPNAMES = {{'Growth','Static'},{'Additive (no \gamma)','Multiplicative'},...
    {'Additive','Additive (no \gamma)'},{'Additive','Multiplicative'},{'Additive','Additive (no \gamma)'}};

MdlTypesInd = [1 2 2 3 3 3 3 3 4 4 4 4 4];

for i = 1:NCompareFigs

%     if ismember(i,[1 5])
% figure('Position',[10 364 1875 566])        
%     else
% figure('Position',[233 364 1388 566])
%     end
figure('Position',[233 364 1388 566])
ax1 = axes('Position',[0.0891    0.2915    0.9015    0.5566]);
[ordered_mdl_data{i},V,topo_mdl_order{i},issig,pval] = PlotMdlResults([rand200_topo_mdls{c1(i,1),c1(i,2)}.Fcv; rand200_topo_mdls{c2(i,1),c2(i,2)}.Fcv],TOPO_MDL_LABELS,'GrpNames',GROUPNAMES{i},...
'DataLabel','\it{F_{CV}}','SigLvl',NStatisticalComparisons,'MdlTypesInd',MdlTypesInd,'MdlTypesNames',MDLTYPES);

a = annotation('textbox',[0.0055    0.9468    0.0588    0.0885],'String',FigureAnnot{i},'EdgeColor','none','FontSize',48);

if ~ismember(i,[1 5])
print([FIGURE_LOCATION,'/Figure',FigureNames{i},FigureAnnot{i},'.png'],'-dpng','-r300')
else
   exportgraphics(gcf,[FIGURE_LOCATION,'/Figure',FigureNames{i},FigureAnnot{i},'.png'],'resolution',300) 
end
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
[phys_mdl_out,~,order_phys_mdls,~,~,lgd1,lgd2] = PlotMdlResults([GROWTH; STATIC],LABELS,'GrpNames',{'Growth','Static'},...
'DataLabel','\it{F_{CV}}','SigLvl',phys_nStatisticalComparisons,'MdlTypesInd',[1 3 3 4 4 5 5 6 6 2],'MdlTypesNames',PHYS_MDLTYPES);

lgd1.Position(2) = 0.91;
lgd2.Position(2) = 0.8399;
ax1.Position(2) = 0.3717;
ax1.XLabel.Position = [10.5000 -0.17 1];
%print([FIGURE_LOCATION,'/Figure4.png'],'-dpng','-r300')
exportgraphics(gcf,[FIGURE_LOCATION,'/Figure4.png'],'resolution',300)

close all

%% Figure 5
% This will take a very very long time to run.
% 
MdlColors = [152,78,163;...
55,126,184;...
228,26,28;...
77,175,74;...
255,127,0;...
66,224,245;...
166,86,40;...
77,190,238]./255;

load('random200_data4physmdl.mat')
for i = 1:length(adjs)

    A = double(adjs{i}>0);
    EmpData{1}(i,:) = sum(A); 
    EmpData{3}(i,:) = betweenness_bin(A);
    EmpData{2}(i,:) = clustering_coef_bu(A)';
    EmpData{4}(i,:) = (sum(A_dist.*A))./EmpData{1}(i,:); 

end

rand200_phys_growth_topo = load('random200_FCV_PhysMdls_Growth_1_output_Topography.mat');

MDL = [4 5 10 1];

Model = rand200_phys_growth_topo.Topography{4};

CorrDATA = rand200_phys_mdls{2}.TopographyCorrs(MDL);
CorrDATA{2} = rand200_phys_mdls{1}.TopographyCorrs{5};

SurfaceData.vertices = lh_inflated_verts;
SurfaceData.faces = lh_faces;
SurfaceData.parc = lh_rand200;

load('Missing_rois.mat')
SurfaceData.MissingROI = rand200_missing;

plotSpatialEmbeddings(EmpData,Model,CorrDATA,SurfaceData,LABELS(MDL),MdlColors([4 4 2 1],:),[FIGURE_LOCATION,'/Figure5']);

close all

%% Figure S2

rand200_data = load('random200_data4topomdl.mat');
add3_optim_mdls = load('random200_TopoMdls_add3_Growth_0_output.mat', 'OptimMdl');
add2_optim_mdls = load('random200_TopoMdls_add2_Growth_0_output.mat', 'OptimMdl');
mult2_optim_mdls = load('random200_TopoMdls_mult2_Growth_0_output.mat', 'OptimMdl');

NETS{1} = mult2_optim_mdls.OptimMdl{3}.min_maxKS.repeats.nets;
NETS{2} = add2_optim_mdls.OptimMdl{3}.min_maxKS.repeats.nets;
NETS{3} = add3_optim_mdls.OptimMdl{3}.min_maxKS.repeats.nets;

for j = 1:3
    IND = 1;
for i = 1:length(NETS{j})
    for k = 1:100
    a = zeros(100);
    a(NETS{j}{k}{i}) = 1;
    a = a + a';
    NETS_FULL{j}{IND} = a;
    IND = IND + 1;
    end
end
end

figure('Position',[680   320   783   658])
CompareEdgeLengtheCDF(rand200_data.A_dist,NETS_FULL{1},NETS_FULL{2},NETS_FULL{3},rand200_data.adjs)

exportgraphics(gcf,[FIGURE_LOCATION,'/FigureS2.png'],'resolution',300)

close all

%% Figure S4

VisualiseKSstats(rand200_topo_mdls{1,1}.KS,TOPO_MDL_LABELS,3,5)
a = annotation('textbox',[0.0018    0.9312    0.0588    0.0885],'String','A','EdgeColor','none','FontSize',48);

print([FIGURE_LOCATION,'/FigureS4A.png'],'-dpng','-r300')

VisualiseKSstats(rand200_topo_mdls{1,2}.KS,TOPO_MDL_LABELS,3,5)
a = annotation('textbox',[0.0018    0.9312    0.0588    0.0885],'String','B','EdgeColor','none','FontSize',48);

print([FIGURE_LOCATION,'/FigureS4B.png'],'-dpng','-r300')

close all

%% Figure S5

% For eta parameters, times by -1, because the paper specifies eta is
% always <0 but the code can accept positive values. To "avoid" confusion,
% in the code eta is always explicitly defined as a negative value, we just
% report it as a positive value in the paper because in the equation, eta
% is always expressed as being converted to a negative value.

for i = 1:13
    
   alphavals_static_add3{i} = rand200_topo_mdls{1,1}.P{i}(:,3);
   alphavals_growth_add3{i} = rand200_topo_mdls{1,2}.P{i}(:,3);
   
   alphavals_static_add2{i} = rand200_topo_mdls{2,1}.P{i}(:,3);
   alphavals_growth_add2{i} = rand200_topo_mdls{2,2}.P{i}(:,3);
 
   etavals_static_add3{i} = rand200_topo_mdls{1,1}.P{i}(:,1)*-1;
   etavals_growth_add3{i} = rand200_topo_mdls{1,2}.P{i}(:,1)*-1; 
   
   gamvals_static_add3{i} = rand200_topo_mdls{1,1}.P{i}(:,2);
   gamvals_growth_add3{i} = rand200_topo_mdls{1,2}.P{i}(:,2);
 
   etavals_static_add2{i} = rand200_topo_mdls{2,1}.P{i}(:,1)*-1;
   etavals_growth_add2{i} = rand200_topo_mdls{2,2}.P{i}(:,1)*-1; 
   
   gamvals_static_mult2{i} = rand200_topo_mdls{3,1}.P{i}(:,2);
%   gamvals_growth_mult2{i} = rand200_topo_mdls{3,2}.P{i}(:,2);
 
   etavals_static_mult2{i} = rand200_topo_mdls{3,1}.P{i}(:,1)*-1;
%   etavals_growth_mult2{i} = rand200_topo_mdls{3,2}.P{i}(:,1)*-1;   
   
end
figure('Position',[233 364 1388 566])
ax1 = axes('Position',[0.0891    0.2915    0.9015    0.5566]);
PlotMdlResults([etavals_growth_add3; etavals_static_add3],TOPO_MDL_LABELS,'GrpNames',{'Growth','Static'},...
'DataLabel','\eta','SigLvl',0,'MdlTypesInd',MdlTypesInd,'MdlTypesNames',MDLTYPES,'MdlOrder',topo_mdl_order{1});
a = annotation('textbox',[0.0055    0.9468    0.0588    0.0885],'String','A','EdgeColor','none','FontSize',48);

%exportgraphics(gcf,[FIGURE_LOCATION,'/FigureS4A.png'],'resolution',300)
print([FIGURE_LOCATION,'/FigureS5A.png'],'-dpng','-r300')

gamvals_static_add3{1} = [];
gamvals_growth_add3{1} = [];

figure('Position',[233 364 1388 566])
ax1 = axes('Position',[0.0891    0.2915    0.9015    0.5566]);

PlotMdlResults([gamvals_growth_add3; gamvals_static_add3],TOPO_MDL_LABELS,'GrpNames',{'Growth','Static'},...
'DataLabel','\gamma','SigLvl',0,'MdlTypesInd',MdlTypesInd,'MdlTypesNames',MDLTYPES,'MdlOrder',topo_mdl_order{1});
a = annotation('textbox',[0.0055    0.9468    0.0588    0.0885],'String','B','EdgeColor','none','FontSize',48);

%exportgraphics(gcf,[FIGURE_LOCATION,'/FigureS4B.png'],'resolution',300)

print([FIGURE_LOCATION,'/FigureS5B.png'],'-dpng','-r300')

alphavals_static_add3{1} = [];
alphavals_growth_add3{1} = [];

figure('Position',[233 364 1388 566])
ax1 = axes('Position',[0.0891    0.2915    0.9015    0.5566]);

PlotMdlResults([alphavals_growth_add3; alphavals_static_add3],TOPO_MDL_LABELS,'GrpNames',{'Growth','Static'},...
'DataLabel','\alpha','SigLvl',0,'MdlTypesInd',MdlTypesInd,'MdlTypesNames',MDLTYPES,'MdlOrder',topo_mdl_order{1});

a = annotation('textbox',[0.0055    0.9468    0.0588    0.0885],'String','C','EdgeColor','none','FontSize',48);

%exportgraphics(gcf,[FIGURE_LOCATION,'/FigureS4C.png'],'resolution',300)

print([FIGURE_LOCATION,'/FigureS5C.png'],'-dpng','-r300')

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
[~,~,~,~,~,lgd1,lgd2]= PlotMdlResults([GROWTH; STATIC],LABELS,'GrpNames',{'Growth','Static'},...
'DataLabel','\it{F_{CV}}','SigLvl',phys_nStatisticalComparisons,'MdlTypesInd',[1 3 3 4 4 5 5 6 6 2],'MdlTypesNames',PHYS_MDLTYPES);

lgd1.Position(2) = 0.91;
lgd2.Position(2) = 0.8399;
ax1.Position(2) = 0.3717;
ax1.XLabel.Position = [10.5000 -.065 1];

exportgraphics(gcf,[FIGURE_LOCATION,'/FigureS7.png'],'resolution',300)

close all

%% Figure S8

rand200_phys_optim{1} = load('random200_PhysMdls_Growth_0_output.mat','DegCorr');
rand200_phys_optim{2} = load('random200_PhysMdls_Growth_1_output.mat','DegCorr');
figure('Position',[233 364 1388 605])
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

[CorrDataOrdered,~,~,~,~,lgd1,lgd2] = PlotMdlResults([GROWTH_DEGCORR; STATIC_DEGCORR],LABELS,'GrpNames',{'Growth','Static'},...
'DataLabel',{'Spearman correlation','with empirical degree'},'SigLvl',0,'MdlTypesInd',[1 3 3 4 4 5 5 6 6 2],'MdlTypesNames',PHYS_MDLTYPES,'MdlOrder',order_phys_mdls);

lgd1.Position(2) = 0.91;
lgd2.Position(2) = 0.8399;
ax1.Position(2) = 0.3717;
ax1.XLabel.Position = [10.5000 -1.714 1];

exportgraphics(gcf,[FIGURE_LOCATION,'/FigureS8.png'],'resolution',300)

close all

%% Figures not in paper
% These were just useful things for reference:

% [h,p] = ComputeSigDiff([rand200_topo_mdls{1,1}.Fcv rand200_topo_mdls{1,2}.Fcv],.05,325,1);
% 
% imagesc(h)
% xticks(1:26)
% yticks(1:26)
% 
% for i = 1:13
%    KS_PLOT_TOPO_LABELS{i} = ['Static ',TOPO_MDL_LABELS{i}];
% end
% for i = 14:26
%    KS_PLOT_TOPO_LABELS{i} = ['Growth ',TOPO_MDL_LABELS{i-13}];
% end
% xticklabels(KS_PLOT_TOPO_LABELS)
% yticklabels(KS_PLOT_TOPO_LABELS)
% xtickangle(45)
% 
% [h,p] = ComputeSigDiff([rand200_topo_mdls{2,1}.Fcv rand200_topo_mdls{2,2}.Fcv],.05,325,1);
% 
% imagesc(h)
% xticks(1:26)
% yticks(1:26)
% 
% for i = 1:13
%    KS_PLOT_TOPO_LABELS{i} = ['Static ',TOPO_MDL_LABELS{i}];
% end
% for i = 14:26
%    KS_PLOT_TOPO_LABELS{i} = ['Growth ',TOPO_MDL_LABELS{i-13}];
% end
% xticklabels(KS_PLOT_TOPO_LABELS)
% yticklabels(KS_PLOT_TOPO_LABELS)
% xtickangle(45)
