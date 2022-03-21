function Compile_PhysMdls_FCV(DATALOC,SAVELOC)

% This compiles the FCV/crossvalidated data for the physiological models
% Inputs:
%
%                   DATALOC = location of the Crossvalidated directory
%
%                   SAVELOC = where to save the outputs

parcnames = {'random200','Schaefer200'};

for parc = 1:length(parcnames)
    for GROWTH = 0:1

    fileformat = [DATALOC,'/Crossvalidated/CrossValidate_',parcnames{parc},'_PhysMdls_mdl_#_Growth_',num2str(GROWTH),'_iter_$.mat'];

    data = CompileCVOutputs(fileformat,1:10,20,100,0);
    Topography = data.Topography;
    data = rmfield(data,'Topography');

    save([SAVELOC,'/',parcnames{parc},'_FCV_PhysMdls_Growth_',num2str(GROWTH),'_output.mat'],'-struct','data','-v7.3');
    save([SAVELOC,'/',parcnames{parc},'_FCV_PhysMdls_Growth_',num2str(GROWTH),'_output_Topography.mat'],'Topography','-v7.3');

    end
end