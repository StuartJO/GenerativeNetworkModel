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

   load(filename)
   N = length(E);
   for j = 1:N
     tic
    [Mdlouts.MinMdlFit{j}(IND),idx] = min(E{j});
    [Mdlouts.MaxCorrFit{j}(IND),idx2] = max(C{j});
    Mdlouts.MdlFits{j}{IND} = E{j};
    Mdlouts.CorrFits{j}{IND} = C{j};
    Mdlouts.KSvals{j}{IND} = K{j};
    Mdlouts.Params{j}{IND} = P{j};
    %alphavals(IND,j) = P{j}(idx,3);
    Mdlouts.MinMdlFitParams{j}(IND,:) = P{j}(idx,:);
    Mdlouts.MaxCorrFitParams{j}(IND,:) = P{j}(idx2,:);
    %alphavals_growth(IND,j) = P{j}(idx,2);
    %etavals_growth(IND,j) = P{j}(idx,1);
    Mdlouts.MinMdlFit_Corr{j}(IND) = C{j}(idx);

    Mdlouts.MinMdlFitParams_MdlFits{j}{IND} = Ebest{j};
    Mdlouts.MinMdlFitParams_Corr{j}{IND} = Cbest{j};
    
	if savenets
		Mdlouts.MinMdlFitParams_Net{j}{IND} = B;
	end

    Asim = zeros(Nodes);
    Asim(b{j}{idx}) = 1;
    Asim = Asim + Asim';
    Mdlouts.MinMdlFitNet{j}{IND} = Asim;
    
    
    Asim = zeros(Nodes);
    Asim(b{j}{idx2}) = 1;
    Asim = Asim + Asim';
    Mdlouts.MaxCorrFitNet{j}{IND} =  Asim;
    
   timesecs = toc;
   disp(['Completed subject ',num2str(i),' model ',num2str(j),' in ',num2str(timesecs),' seconds'])
   end
    
   IND = IND + 1;
   
    else
        
       Mdlouts.missing = [Mdlouts.missing i];
    end
      
end