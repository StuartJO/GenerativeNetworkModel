function Compile_TopoMdls_FCV(DATALOC,SAVELOC)

% This compiles the FCV/crossvalidated data for the topological models
% Inputs:
%
%                   DATALOC = location of the Crossvalidated directory
%
%                   SAVELOC = where to save the outputs

parcs = {'random200'};
types = {'mult2','add2','add3'};
for parc = 1:length(parcs)
    for type = 1:3
        if type == 1
            GROWTHVALS = 0;
        else
            GROWTHVALS = 0:1;
        end
        for GROWTH = GROWTHVALS

        fileformat = [DATALOC,'/Crossvalidated/CrossValidate_',parcs{parc},'_TopoMdls_',types{type},'_mdl_#_Growth_',num2str(GROWTH),'_iter_$.mat'];

        data = CompileCVOutputs(fileformat,1:13,20,100,0);
        Topography = data.Topography;
        data = rmfield(data,'Topography');

        save([SAVELOC,'/',parcs{parc},'_FCV_TopoMdls_',types{type},'_Growth_',num2str(GROWTH),'_output.mat'],'-struct','data','-v7.3');
        save([SAVELOC,'/',parcs{parc},'_FCV_TopoMdls_',types{type},'_Growth_',num2str(GROWTH),'_output_Topography.mat'],'Topography','-v7.3');
        end
    end
end