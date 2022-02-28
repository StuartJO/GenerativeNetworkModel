function Mdlouts = CompileGenMdl(fileformat,Nodes,savenets)

if nargin < 3
	savenets = 0;
end

Mdlouts.missing = [];

filesubind = strfind(fileformat,'$');
file_start = fileformat(1:filesubind-1);
file_end = fileformat(filesubind+1:end);

IND = 1;
for i = 1:100
filename = [file_start,num2str(i),file_end];
if exist(filename) == 2

    filedata = load(filename);
    Nmdls = length(filedata.E);
        for j = 1:Nmdls
            
        tic
        
        [Mdlouts.min_maxKS{j}(IND),min_maxKS_ind] = min(filedata.E{j});
        [Mdlouts.max_DegCorr{j}(IND),max_DegCorr_ind] = max(filedata.C{j});
        Mdlouts.maxKS{j}{IND} = filedata.E{j};
        Mdlouts.DegCorr{j}{IND} = filedata.C{j};
        Mdlouts.KS{j}{IND} = filedata.K{j};
        Mdlouts.P{j}{IND} = filedata.P{j};

        Mdlouts.OptimMdls.min_maxKS.min_maxKS{j}(IND,:) = min(filedata.E{j});
        Mdlouts.OptimMdls.max_DegCorr.max_DegCorr{j}(IND,:) = max(filedata.C{j}); 

        Mdlouts.OptimMdls.min_maxKS.P{j}(IND,:) = filedata.P{j}(min_maxKS_ind,:);
        Mdlouts.OptimMdls.max_DegCorr.P{j}(IND,:) = filedata.P{j}(max_DegCorr_ind,:);
        Mdlouts.OptimMdls.min_maxKS.DegCorr{j}(IND) = filedata.C{j}(min_maxKS_ind);
        Mdlouts.OptimMdls.max_DegCorr.maxKS{j}(IND) = filedata.E{j}(max_DegCorr_ind);

        Mdlouts.OptimMdls.min_maxKS.repeats.maxKS{j}{IND} = filedata.Ebest{j};
        Mdlouts.OptimMdls.min_maxKS.repeats.DegCorr{j}{IND} = filedata.Cbest{j};

            if savenets == 1

            Asim = zeros(Nodes);
            Asim(filedata.b{j}{min_maxKS_ind}) = 1;
            Asim = Asim + Asim';
            Mdlouts.OptimMdls.min_maxKS.repeats.nets{j}{IND} = Asim;

            Asim = zeros(Nodes);
            Asim(filedata.b{j}{max_DegCorr_ind}) = 1;
            Asim = Asim + Asim';
            Mdlouts.OptimMdls.max_DegCorr.repeats.nets{j}{IND} =  Asim;

            elseif savenets == 2
                
            Asim = zeros(Nodes);
            Asim(filedata.b{j}{min_maxKS_ind}) = 1;
            Asim = Asim + Asim';
            Mdlouts.OptimMdls.min_maxKS.repeats.nets{j}{IND} = Asim;

            Asim = zeros(Nodes);
            Asim(filedata.b{j}{max_DegCorr_ind}) = 1;
            Asim = Asim + Asim';
            Mdlouts.OptimMdls.max_DegCorr.repeats.nets{j}{IND} =  Asim;
            Mdlouts.nets_edgeind{j}{IND} = filedata.b;

            end

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