function Compile_PhysMdls(DATALOC,SAVELOC)

% This compiles the optimisation data for the physiological models
% Inputs:
%
%                   DATALOC = location of the optimisation directory
%
%                   SAVELOC = where to save the outputs

parcnames = {'random200','Schaefer200'};
for i = 1:2
    for GROWTH = 0:1
        fileformat = [DATALOC,'/Optimisation/',parcnames{i},'_PhysMdls_Sub_#_add3_Growth_',num2str(GROWTH),'.mat'];
        Mdlouts = CompileGenMdlOutputs(fileformat,1:100,1);
              
        save([SAVELOC,parcnames{i},'_PhysMdls_Growth_',num2str(GROWTH),'_output.mat'],'-struct','Mdlouts','-v7.3');  
    end
end