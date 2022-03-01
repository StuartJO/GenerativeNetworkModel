function CV = Compile_FCV(fileformat,mdls,Iters,nsub,SaveNets)

% This function allows you to compile CV results into a single output.
%
% Inputs:
%
% fileformat = a string giving the format of the file to load in.
% Because the output of CV scripts is done for different models and
% different iterations, you can use '#' to signify a model number in the
% file format, and '$' to indicate the iteration number. For example:
%
% Note that this script assumes each model is its own file. It is optional
% as to whether many files of different iterations have been made (i.e., if
% all iterations were saved to the one file already, this script will
% handle that, just don't include '$' in the fileformat). Also this script
% assumes the model number comes before the iteration number
%
% mdls = an array indicating which number models to include, and the order.
% mdls = [1 5 4] woudl include the models numbers 1, 5 and 4 in that order
%
% Iters = a scaler or 1*2 array. The scaler indicates the number of
% iterations in each file if '$' is included in the fileformat, or the
% number of iteration files otherwise. If an array, the first value
% indiates the number of iteration files and the second the number of 
% iterations in each file
%
% nsub = the number of subjects for whom cross validation was performed.
%
% SaveNets = set to 1 to save all the networks produced by crossvalidation
% on top of all the statistic (0 by default).
%
% Output:

if nargin < 5
    SaveNets = 0;
end

% Checks if there are multiple iteration files
if ~contains(fileformat,'$')
    % Multiple iteration files do not exist
    NoIterFiles = 1;

    % divide the fileformat string into sections based on the position of
    % the '#' character
    filesubind1 = strfind(fileformat,'#');
    file_start = fileformat(1:filesubind1-1);

    file_end = fileformat(filesubind1+1:end);
    
    nperIterFile = Iters;
    IterFiles = 1;
else
    % Multiple iteration files do exist
    NoIterFiles = 0;

    % divide the fileformat string into sections based on the position of
    % the '#' and '$' characters
    filesubind1 = strfind(fileformat,'#');
    file_start = fileformat(1:filesubind1-1);

    filesubind2 = strfind(fileformat,'$');
    
    if filesubind1 > filesubind2
       error('model number needs to come before iteration number in the filename!') 
    end
    
    file_middle = fileformat(filesubind1+1:filesubind2-1);
    file_end = fileformat(filesubind2+1:end);
    
    if length(Iters) == 1
    nperIterFile = 1;
    IterFiles = 1:Iters;
    else
    nperIterFile = Iters(2);
    IterFiles = 1:Iters(1);        
    end
end

nIters = length(IterFiles)*nperIterFile;

totalFcvDataPoints = ((nsub*nsub)-nsub)*nIters;

Mdsl2Compile = length(mdls);

for j = 1:Mdsl2Compile
    tic
    PointsPerIter = ((nsub*nsub)-nsub)*nperIterFile;
    data = zeros(nsub,nsub,nIters);
    if SaveNets
        nets = cell(nsub,nsub);
    end
    CV.KS{j} = zeros(totalFcvDataPoints,4);
    CV.TopographyCorrs{j} = zeros(totalFcvDataPoints,4);
    CV.DataID{j} = zeros(totalFcvDataPoints,3);
    Topography = cell(1,4);
    
    ind = 1:nperIterFile;
     
    for itr = 1:length(IterFiles)
        if NoIterFiles
           CVindata = load([file_start,num2str(mdls(j)),file_end]);
        else
           CVindata = load([file_start,num2str(mdls(j)),file_middle,num2str(itr),file_end]) ;   
        end
        
        NNodes = CVindata.Input.NNodes;
        
    for t = 1:4
          Topography{t}{itr} = zeros(PointsPerIter,NNodes);
    end          
   topo_ind = 1:nperIterFile;
       % sub_n indicates that index corresponds to the subject network
       for sub_n = 1:nsub
           %sub_n indicates that index corresponds to the subject parameter
            for sub_P = 1:nsub
            if sub_n == sub_P
                % CV should not include the subjects own parameters applied to
                % their own network. We set the respective indicies to nan
                data(sub_n,sub_P,(1:nperIterFile)+(nperIterFile*(itr-1))) = nan(1,nperIterFile);       
            else
                data(sub_n,sub_P,(1:nperIterFile)+(nperIterFile*(itr-1))) = CVindata.KSmax{sub_n,sub_P};
                
 % You might wonder why I don't directly put the topography data direct
 % into the output structure, like I do with all the other outputs. The
 % reason is simple/dumb. For whatever reason, when I try do it
 % directly, this function slows down by an order of magnitude. Even if
 % I initalise Fcv.Topography, it still does not work. However, if I define
 % a variable (in this case "Topography"), write to that, and then
 % put that into the output structure, well that just works fine. This is 
 % very much a case of "If it is stupid but it works...."                
                
                for t = 1:4
                    Topography{t}{itr}(topo_ind,:) = CVindata.Topography{sub_n,sub_P}{t}(1:nperIterFile,:);
                end 
                
                % Save all the KS and correlation data into a huge array
                CV.KS{j}(ind,:) = CVindata.KS{sub_n,sub_P}(1:nperIterFile,:);
                CV.TopographyCorrs{j}(ind,:) = CVindata.TopographyCorrs{sub_n,sub_P}(1:nperIterFile,:);
                % CV.DataID gives a "barcode" indicating which network,
                % which set of parameters and for which iteration that
                % respective data belongs to
                CV.DataID{j}(ind,:) = [sub_n sub_P itr];

                ind = ind + nperIterFile;
                topo_ind = topo_ind+nperIterFile;
            end
            
            if SaveNets
                nets{sub_n,sub_P}((1:nperIterFile)+(nperIterFile*(itr-1))) = CVindata.Net{sub_n,sub_P};
            end
            
           end

       end
   
    end
     CV.Inputs{j} = CVindata.Input;
% "...it ain't stupid"     
     topo_ind2 = 1:PointsPerIter;
     for itr = 1:length(IterFiles)
         for t = 1:4
            CV.Topography{j}{t}(topo_ind2,:) = Topography{t}{itr};
         end
        topo_ind2 = topo_ind2+PointsPerIter;
     end
    clear Topography

        timesec = toc;
        disp(['Finished ',num2str(j),' in ',num2str(timesec),' seconds'])
        
    if SaveNets
    	CV.nets{j} = nets;
    end
    % Take the mean over all iterations, and then take the mean over all
    % the parameters applied to a given subject
   CV.Fcv{j} = nanmean(nanmean(data,3),2)'; 
   CV.P{j} = CVindata.P;
end