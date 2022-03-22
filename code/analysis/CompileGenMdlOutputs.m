function Mdlouts = CompileGenMdlOutputs(fileformat,subs,savenets,MltpFileINDS)

% This function compiles outputs (i.e., different .mat files) for a given
% set of models and saves it to one convenient structure.
%
% Inputs:
%
% fileformat = a string giving the format of the file to load in.
% You can use '#' to signify each subjects ID, and optionally '$' to 
% indicate the model number. For example:
%
%   fileformat = 'random200_TopoMdls_Sub_#_add3_Growth_0.mat'
%
%   This will loop over the values in the variable "subs", where '#' will
%   correspond to each value.
%
%   fileformat = 'random200_TopoMdls_Sub_#_add3_Growth_0_mdl_$.mat'
%
%   As before, but additonally this will loop over the values in the
%   variable "MltpFileINDS", where '$' will correspond to each value
%
%
% The code implicitly assumes each file has the same number of models in
% it. Also this script assumes the subject ID comes before the model number
%
% subs = array of subject IDs to loop over to load in
%
% savenets = set to 1 to save only the models generated for repeats of the
% best performing models, set to 2 to also save all networks produced
% during the optimisation
%
% MltpFileINDS = a list of indices to loop over the '$' value in 
% fileformat (only required if '$' is specified)
% 
% Outputs:
% 
% Mdlouts = a structure containing the following fields:
%
%   min_maxKS{mdlIND}(subIND) = the minimum maxKS for model 'mdlIND' and 
%       subject 'subIND' across all networks produced during optimisation
%       
%   max_DegCorr{mdlIND}(subIND) = the maximum correlation with empirical 
%       degree for model 'mdlIND' and subject 'subIND' across all networks 
%       produced during optimisation
%
%   maxKS{mdlIND}{subIND} = the maxKS for model 'mdlIND' and subject 'subIND' 
%       for all networks produced during optimisation
%
%   DegCorr{mdlIND}{subIND} = the correlation with empirical 
%       degree for model 'mdlIND' and subject 'subIND' for all networks 
%       produced during optimisation
%
%   KS{mdlIND}{subIND} = the KS values for  model 'mdlIND' and subject 
%       'subIND' for all networks produced during optimisation
%       
%   P{mdlIND}{subIND} = the parameters for model 'mdlIND' and subject 
%       'subIND' (N*5 matrix, where N is the number of networks produced 
%       during optimisation)  
% 
%   OptimMdl{mdlIND}.min_maxKS.min_maxKS(subIND,:) = the minimum maxKS for 
%       model 'mdlIND' and subject 'subIND' across all networks produced
%       during optimisation
%
%   OptimMdl{mdlIND}.min_maxKS.P(subIND,:) = for model 'mdlIND' and 
%       subject 'subIND', the parameters which produced the network with
%       the smallest maxKS during optimisation
%
%   OptimMdl{mdlIND}.min_maxKS.DegCorr(subIND) = for the network with the
%       minimum maxKS produced during optimisation for model 'mdlIND' and 
%       subject 'subIND', its correlation with empirical degree
%       
%   OptimMdl{mdlIND}.min_maxKS.repeats.maxKS{subIND} = for all networks
%       produced using the parameters which produced the network with the
%       minimum maxKS for model 'mdlIND' and subject 'subIND', the maxKS of
%       those
%
%   OptimMdl{mdlIND}.min_maxKS.repeats.DegCorr{subIND} = for all networks
%       produced using the parameters which produced the network with the
%       minimum maxKS for model 'mdlIND' and subject 'subIND', the
%       correlation with empirical degree for those
%
%   OptimMdl{mdlIND}.min_maxKS.net{subIND} = cell containing the network
%       the produced the lowestest maxKS for model 'mdlIND' and subject 'subIND'
%
%   BestMdls{mdlIND}.max_DegCorr.max_DegCorr(subIND,:) = the maximum 
%       correlation with empirical degree for model 'mdlIND' and subject
%        'subIND' across all networks produced during optimisation
% 
%   BestMdls{mdlIND}.max_DegCorr.P(subIND,:) = for model 'mdlIND' and 
%       subject 'subIND', the parameters which produced the network with
%       the largest correlation with empirical degree
%
%   BestMdls{mdlIND}.max_DegCorr.maxKS(subIND) = for the network with the
%       largest correlation with empirical degree produced during 
%       optimisation for model 'mdlIND' and subject 'subIND', its maxKS
%       value
% 
%   BestMdls{mdlIND}.max_DegCorr.repeats.maxKS{subIND} = for all networks
%       produced using the parameters which produced the network with the
%       largest correlation with empirical degree for model 'mdlIND' and
%       subject 'subIND', the maxKS for those
%
%   BestMdls{mdlIND}.max_DegCorr.repeats.DegCorr{subIND} = for all networks
%       produced using the parameters which produced the network with the
%       largest correlation with empirical degree for model 'mdlIND' and
%       subject 'subIND', the correlation with empirical degre for those
%
%   BestMdls{mdlIND}.max_DegCorr.net{subIND} = cell containing the network
%       the produced the highest correlation with emprical degree for model
%       'mdlIND' and subject 'subIND'
%
%   nets{mdlIND}{subIND} = all the networks produced during optimisation
%       for model 'mdlIND' and subject 'subIND'
%
%   Input{mdlIND} = the initial input settings. NOTE: this assumes the saem
%   settings for all subjects

if nargin < 2
	savenets = 0;
end

Mdlouts.missing = [];

% Checks if there are multiple model files for a single subject
if ~contains(fileformat,'$')
    % Multiple model files do not exist
    MltpMdlFilesExist = 0;

    % divide the fileformat string into sections based on the position of
    % the '#' character
    filesubind1 = strfind(fileformat,'#');
    file_start = fileformat(1:filesubind1-1);

    file_end = fileformat(filesubind1+1:end);
    
    MltpFileINDS = 1;
    
else
    % Multiple model files for a single subject do exist
    MltpMdlFilesExist = 1;

    % divide the fileformat string into sections based on the position of
    % the '#' and '$' characters
    filesubind1 = strfind(fileformat,'#');
    file_start = fileformat(1:filesubind1-1);

    filesubind2 = strfind(fileformat,'$');
    
    if filesubind1 > filesubind2
       error('subject number needs to come before model number in the filename!') 
    end
    
    file_middle = fileformat(filesubind1+1:filesubind2-1);
    file_end = fileformat(filesubind2+1:end);
    
end

% filesubind = strfind(fileformat,'#');
% file_start = fileformat(1:filesubind-1);
% file_end = fileformat(filesubind+1:end);

subIND = 1;

for sub = subs
    mdlIND = 1;
    for MltpMdlFile = MltpFileINDS
        if MltpMdlFilesExist
           filename = [file_start,num2str(sub),file_middle,num2str(MltpMdlFile),file_end];
        else
           filename = [file_start,num2str(sub),file_end]; 
        end

        if exist(filename) == 2

            filedata = load(filename);
            Nmdls = length(filedata.maxKS);
                for mdl = 1:Nmdls

                tic

                [Mdlouts.min_maxKS{mdlIND}(subIND),min_maxKS_ind] = min(filedata.maxKS{mdl});
                [Mdlouts.max_DegCorr{mdlIND}(subIND),max_DegCorr_ind] = max(filedata.DegCorr{mdl});

                Mdlouts.maxKS{mdlIND}{subIND} = filedata.maxKS{mdl};
                Mdlouts.DegCorr{mdlIND}{subIND} = filedata.DegCorr{mdl};
                Mdlouts.KS{mdlIND}{subIND} = filedata.KS{mdl};
                Mdlouts.P{mdlIND}{subIND} = filedata.P{mdl};

                Mdlouts.OptimMdl{mdlIND}.min_maxKS.min_maxKS(subIND,:) = min(filedata.maxKS{mdl});
                Mdlouts.BestMdls{mdlIND}.max_DegCorr.max_DegCorr(subIND,:) = max(filedata.DegCorr{mdl}); 

                Mdlouts.OptimMdl{mdlIND}.min_maxKS.P(subIND,:) = filedata.P{mdl}(min_maxKS_ind,:);
                Mdlouts.BestMdls{mdlIND}.max_DegCorr.P(subIND,:) = filedata.P{mdl}(max_DegCorr_ind,:);
                Mdlouts.OptimMdl{mdlIND}.min_maxKS.DegCorr(subIND) = filedata.DegCorr{mdl}(min_maxKS_ind);
                Mdlouts.BestMdls{mdlIND}.max_DegCorr.maxKS(subIND) = filedata.maxKS{mdl}(max_DegCorr_ind);

                Mdlouts.OptimMdl{mdlIND}.min_maxKS.repeats.maxKS{subIND} = filedata.optim_maxKS{mdl};
                Mdlouts.OptimMdl{mdlIND}.min_maxKS.repeats.DegCorr{subIND} = filedata.optim_DegCorr{mdl};

                Mdlouts.BestMdls{mdlIND}.max_DegCorr.repeats.maxKS{subIND} = filedata.bestDegCorr_maxKS{mdl};
                Mdlouts.BestMdls{mdlIND}.max_DegCorr.repeats.DegCorr{subIND} = filedata.bestDegCorr_DegCorr{mdl};

                Mdlouts.OptimMdl{mdlIND}.min_maxKS.net{subIND} = filedata.b{mdl}{min_maxKS_ind};
                Mdlouts.BestMdls{mdlIND}.max_DegCorr.net{subIND} = filedata.b{mdl}{max_DegCorr_ind};

                    if savenets == 1

                    Mdlouts.OptimMdl{mdlIND}.min_maxKS.repeats.nets{subIND} = filedata.optim_b{mdl};    
                    Mdlouts.BestMdls{mdlIND}.max_DegCorr.repeats.nets{subIND} = filedata.bestDegCorr_b{mdl}; 

                    elseif savenets == 2

                    Mdlouts.OptimMdl{mdlIND}.min_maxKS.repeats.nets{subIND} = filedata.optim_b{mdl};    
                    Mdlouts.BestMdls{mdlIND}.max_DegCorr.repeats.nets{subIND} = filedata.bestDegCorr_b{mdl};                

                    Mdlouts.nets{mdlIND}{subIND} = filedata.b{mdl};

                    end
                    % I assume that across subjects, the 'Input' variable is the
                    % same for a given model.
                Mdlouts.Inputs{mdlIND} = filedata.Input{mdl}; 
                timesecs = toc;
                disp(['Completed subject ',num2str(sub),' model ',num2str(mdlIND),' in ',num2str(timesecs),' seconds'])

                mdlIND = mdlIND + 1;
                end

            subIND = subIND + 1;

        else
            % Sometimes I would run the script and not all models had finished
            % running so there would be no output. This just allowed me to 
            % identify any such subjects 
            Mdlouts.missing = [Mdlouts.missing sub];
            disp(['Subject ',num2str(sub),' data not found!'])
        end

    end

end