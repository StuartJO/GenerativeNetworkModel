function [Fcv,NETS] = Compile_FCV(fileformat,mdls,Iters,nsub,N)

if ~contains(fileformat,'$')
    NoIterFiles = 1;

filesubind1 = strfind(fileformat,'#');
file_start = fileformat(1:filesubind1-1);

file_end = fileformat(filesubind1+1:end);
nperIter = Iters;
iters2use = 1;
else
    NoIterFiles = 0;

filesubind1 = strfind(fileformat,'#');
file_start = fileformat(1:filesubind1-1);

filesubind2 = strfind(fileformat,'$');
file_middle = fileformat(filesubind1+1:filesubind2-1);
file_end = fileformat(filesubind2+1:end);
nperIter = 1;
iters2use = 1:Iters;
end



addpath /projects/kg98/stuarto/BCT

%nsub = 100;

ntotal = length(iters2use)*nperIter;



totalDataPoints = ((nsub*nsub)-nsub)*Iters;

Mdsl2Compile = length(mdls);

for j = 1:Mdsl2Compile
%     for k = 1:4
% Fcv.Topography{j}{k} = zeros(99*99*20,N);
% end
%     Fcv.DataID{j} = zeros(99*99*20,3);
    tic
tic
PointsPerIter = ((nsub*nsub)-nsub)*nperIter;
     data = zeros(nsub,nsub,ntotal);
     nets = cell(nsub,nsub);
    
%      for k = 1:4
%            Fcv.Topography{j}{k} = zeros(totalDataPoints,N);
%      end 
     Fcv.KS{j} = zeros(totalDataPoints,4);
     Fcv.TopographyCorrs{j} = zeros(totalDataPoints,4);
     Fcv.DataID{j} = zeros(totalDataPoints,3);
     ind = 1:nperIter;
     
     for l = 1:length(iters2use)%1:5
         if NoIterFiles
            CVdata = load([file_start,num2str(mdls(j)),file_end]);
         else
            CVdata = load([file_start,num2str(mdls(j)),file_middle,num2str(l),file_end]) ;   
         end
     for k = 1:4
           Topography{k}{l} = zeros(PointsPerIter,N);
     end          
    topo_ind = 1:nperIter;
   % i is the subjects network
   for i = 1:nsub
       %ii is parameter
       for ii = 1:nsub
        if i == ii
                data(i,ii,[1:nperIter]+(nperIter*(l-1))) = nan(1,nperIter);
               
        else
                data(i,ii,[1:nperIter]+(nperIter*(l-1))) = CVdata.KSmax{i,ii};
              for k = 1:4
           %Fcv.Topography{j}{k}(ind,:) = CVdata.Topography{i,ii}{k}(1:nperIter,:);
           Topography{k}{l}(topo_ind,:) = CVdata.Topography{i,ii}{k}(1:nperIter,:);
              end 
                Fcv.KS{j}(ind,:) = CVdata.KS{i,ii}(1:nperIter,:);
                Fcv.TopographyCorrs{j}(ind,:) = CVdata.TopographyCorrs{i,ii}(1:nperIter,:);
                Fcv.DataID{j}(ind,:) = [i ii l];
                
                ind = ind + nperIter;
                topo_ind = topo_ind+nperIter;
        end
        
        nets{i,ii}([1:nperIter]+(nperIter*(l-1))) = CVdata.Net{i,ii};
        
       end
        
   end
   
     end
     
     topo_ind2 = 1:PointsPerIter;
     for l = 1:length(iters2use)
         for k = 1:4
     Fcv.Topography{j}{k}(topo_ind2,:) = Topography{k}{l};
         end
     topo_ind2 = topo_ind2+PointsPerIter;
     end
    clear Topography

%	[minIter,minIterInd] = nanmin(data,[],3);
% 
%    for i = 1:nsub
% 
% 	[~,BEST_PARAM] = nanmin(minIter(i,:));
% 	
% 	[~,BEST_NET_IND] = nanmin(squeeze(data(i,BEST_PARAM,:)));
% 
% 	Fcv.bestnet{j}{i} = nets{i,BEST_PARAM}{BEST_NET_IND};
% 
%    end
%    
%    display(['Finding best networks'])
% 
%    bestNetPerParam = cell(1,nsub);
%   for i = 1:100
%     ITER = 1;
%     for ii = 1:100
%         if ii ~= i
%         	bestNetPerParam{i}{ITER} = nets{i,ii}{minIterInd(i,ii)};
%         	ITER = ITER + 1;
%         end
%     end
%   end





        timesec = toc;
        disp(['Finished ',num2str(j),' in ',num2str(timesec),' seconds'])

   NETS{j} = nets;
   Fcv.SUBDATA{j} = data;
   %Fcv.corr{j} = cdata;
   %Fcv.minmin{j} = nanmin(nanmin(data,[],3),[],2)';
   Fcv.Fcv{j} = nanmean(nanmean(data,3),2)'; 
   %Fcv.meanmin{j} = nanmean(nanmin(data,[],3),2)';
   %Fcv.minmean{j} = nanmin(nanmean(data,3),[],2)'; 
   %Fcv.MinParamNet{j} = bestNetPerParam;
   Fcv.P{j} = CVdata.P;
   
   %Fcv.corr_maxmax{j} = nanmax(nanmax(cdata,[],3),[],2)';
   %Fcv.Fcv_corr{j} = nanmean(nanmean(cdata,3),2)'; 
   %Fcv.corr_meanmax{j} = nanmean(nanmax(cdata,[],3),2)';
   %Fcv.corr_maxmean{j} = nanmax(nanmean(cdata,3),[],2)'; 
end