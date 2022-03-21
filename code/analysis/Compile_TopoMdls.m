function Compile_TopoMdls(DATALOC,SAVELOC)

% This compiles the optimisation data for the topological models
% Inputs:
%
%                   DATALOC = location of the optimisation directory
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

        fileformat = [DATALOC,'/Optimisation/',parcs{parc},'_TopoMdls_Sub_#_',types{type},'_Growth_',num2str(GROWTH),'.mat'];
        
        Mdlouts = CompileGenMdlOutputs(fileformat,1:100,1);
        
        save([SAVELOC,'/',parcs{parc},'_TopoMdls_',types{type},'_Growth_',num2str(GROWTH),'_output.mat'],'-struct','Mdlouts','-v7.3');

        end
    end
end