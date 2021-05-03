function Mdlouts = CompileGenMdl(fileformat,Nodes,D,savenets)

if nargin < 4
	savenets = 0;
end

Mdlouts.missing = [];

%load('random200_100randsubs_for_modelling.mat')

filesubind = strfind(fileformat,'$');
file_start = fileformat(1:filesubind-1);
file_end = fileformat(filesubind+1:end);

IND = 1;
for i = 1:100
    %filename = ['Sub_',num2str(i),'_',filesuffix,'.mat'];
filename = [file_start,num2str(i),file_end];
    if exist(filename) == 2

    %A = double(adjs{i}>0);

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
    
    run = 0;
    if run
    
    deg_vals = zeros(10000,Nodes);
    clu_vals = zeros(10000,Nodes);
    bet_vals = zeros(10000,Nodes);
    dist_vals = zeros(10000,nnz(Asim)/2);
    
    for k = 1:10000
        A = zeros(Nodes);
        A(b{j}{k}) = 1;
        A = A + A';
        deg_vals(k,:) = sum(A,2);
        clu_vals(k,:) = clustering_coef_bu(A);
        bet_vals(k,:) = betweenness_bin(A)';
        dist_vals(k,:) = D(triu(A,1) > 0);
    end
    Mdlouts.MdlNodeDeg{j}{IND} = deg_vals;
    Mdlouts.MdlNodeClu{j}{IND} = clu_vals;
    Mdlouts.MdlNodeBet{j}{IND} = bet_vals;
    Mdlouts.MdlNodeMeanDist{j}{IND} = dist_vals;
    end
    
   timesecs = toc;
   disp(['Completed subject ',num2str(i),' model ',num2str(j),' in ',num2str(timesecs),' seconds'])
   end
    
   IND = IND + 1;
   
    else
        
       Mdlouts.missing = [Mdlouts.missing i];
    end
      
end