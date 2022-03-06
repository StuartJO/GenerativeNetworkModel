function Mdlouts = CompileGenMdlOutputs(fileformat,subs,savenets)

% This function compiles outputs (i.e., different .mat files) for a given
% set of models and saves it to one convenient structure.
%
% Inputs:
%
% fileformat = a string giving the format of the file to load in.
% Use '$' to signify the subject number in the file name. For example:

if nargin < 2
	savenets = 0;
end

Mdlouts.missing = [];

filesubind = strfind(fileformat,'$');
file_start = fileformat(1:filesubind-1);
file_end = fileformat(filesubind+1:end);

IND = 1;
for i = subs
filename = [file_start,num2str(i),file_end];
if exist(filename) == 2

    filedata = load(filename);
    Nmdls = length(filedata.maxKS);
        for j = 1:Nmdls
            
        tic
        
        [Mdlouts.min_maxKS{j}(IND),min_maxKS_ind] = min(filedata.maxKS{j});
        [Mdlouts.max_DegCorr{j}(IND),max_DegCorr_ind] = max(filedata.DegCorr{j});
        
        Mdlouts.maxKS{j}{IND} = filedata.maxKS{j};
        Mdlouts.DegCorr{j}{IND} = filedata.DegCorr{j};
        Mdlouts.KS{j}{IND} = filedata.KS{j};
        Mdlouts.P{j}{IND} = filedata.P{j};

        Mdlouts.OptimMdl{j}.min_maxKS.min_maxKS(IND,:) = min(filedata.maxKS{j});
        Mdlouts.BestMdls{j}.max_DegCorr.max_DegCorr(IND,:) = max(filedata.DegCorr{j}); 

        Mdlouts.OptimMdl{j}.min_maxKS.P(IND,:) = filedata.P{j}(min_maxKS_ind,:);
        Mdlouts.BestMdls{j}.max_DegCorr.P(IND,:) = filedata.P{j}(max_DegCorr_ind,:);
        Mdlouts.OptimMdl{j}.min_maxKS.DegCorr(IND) = filedata.DegCorr{j}(min_maxKS_ind);
        Mdlouts.BestMdls{j}.max_DegCorr.maxKS(IND) = filedata.maxKS{j}(max_DegCorr_ind);

        Mdlouts.OptimMdl{j}.min_maxKS.repeats.maxKS{IND} = filedata.optim_maxKS{j};
        Mdlouts.OptimMdl{j}.min_maxKS.repeats.DegCorr{IND} = filedata.optim_DegCorr{j};

        Mdlouts.BestMdls{j}.max_DegCorr.repeats.maxKS{IND} = filedata.bestDegCorr_maxKS{j};
        Mdlouts.BestMdls{j}.max_DegCorr.repeats.DegCorr{IND} = filedata.bestDegCorr_DegCorr{j};

        Mdlouts.OptimMdl{j}.min_maxKS.net{IND} = filedata.b{j}{min_maxKS_ind};
        Mdlouts.BestMdls{j}.max_DegCorr.net{IND} = filedata.b{j}{max_DegCorr_ind};
                    
            if savenets == 1
            
            Mdlouts.OptimMdl{j}.min_maxKS.repeats.nets{IND} = filedata.optim_b{j};    
            Mdlouts.BestMdls{j}.max_DegCorr.repeats.nets{IND} = filedata.bestDegCorr_b{j}; 

            elseif savenets == 2

            Mdlouts.OptimMdl{j}.min_maxKS.repeats.nets{IND} = filedata.optim_b{j};    
            Mdlouts.BestMdls{j}.max_DegCorr.repeats.nets{IND} = filedata.bestDegCorr_b{j};                

            Mdlouts.nets{j}{IND} = filedata.b;

            end
            % I assume that across subjects, the 'Input' variable is the
            % same for a given model.
        Mdlouts.Inputs{j} = filedata.Input{j}; 
        timesecs = toc;
        disp(['Completed subject ',num2str(i),' model ',num2str(j),' in ',num2str(timesecs),' seconds'])
        end
    
    IND = IND + 1;

else
    % Sometimes I would run the script and not all models had finished
    % running so there would be no output. This just allowed me to 
    % identify any such subjects 
    Mdlouts.missing = [Mdlouts.missing i];
    disp(['Subject ',num2str(i),' data not found!'])
end

end