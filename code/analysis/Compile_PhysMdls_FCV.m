for i = 0:1

fileformat = '/fs02/hf49/Stuart/GrowthModel_newParc/GenerativeNetworkModel/data/Crossvalidated/CrossValidate_random200_PhysMdls_mdl_#_Growth_1_iter_$.mat';

data = Compile_FCV(fileformat,1:10,20,100,0);
Topography = data.Topography;
data = rmfield(data,'Topography');

save(['random200_FCV_PhysMdls_Growth_',num2str(i),'_output.mat'],'-struct','data','-v7.3');
save(['random200_FCV_PhysMdls_Growth_',num2str(i),'_output_Topography.mat'],'Topography','-v7.3');

end