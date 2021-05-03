function [Fcv,NETS] = CompileCrossValidData(fileformat,nIters)

filesubind1 = strfind(fileformat,'#');
file_start = fileformat(1:filesubind1-1);

filesubind2 = strfind(fileformat,'$');
file_middle = fileformat(filesubind1+1:filesubind2-1);
file_end = fileformat(filesubind2+1:end);

addpath /projects/kg98/stuarto/BCT
%mtype = {'sptl','neighbors','matching','clu-avg','clu-min','clu-max','clu-diff','clu-prod','deg-avg','deg-min','deg-max','deg-diff','deg-prod','com'};

nsub = 100;

nperIter = 1;

iters2use = 1:nIters;

ntotal = length(iters2use)*nperIter;

for j = 1:13
    tic
    %modeltype = mtype{j};
    %display(['Loading model ',modeltype])
     data = zeros(nsub,nsub,ntotal);
     cdata = zeros(nsub,nsub,ntotal);  
     nets = cell(nsub,nsub);
     for l = 1:length(iters2use)%1:5
         load([file_start,num2str(j),file_middle,num2str(l),file_end])
   for i = 1:nsub
       for ii = 1:nsub
%         d = mean(EcrossValid{i},2);
%         d(i,:) = NaN;
%         data = [data d];
        if i == ii
                data(i,ii,[1:nperIter]+(nperIter*(l-1))) = nan(1,nperIter);
                cdata(i,ii,[1:nperIter]+(nperIter*(l-1))) = nan(1,nperIter);
        else
                data(i,ii,[1:nperIter]+(nperIter*(l-1))) = EcrossValid{ii}(i,:);
                cdata(i,ii,[1:nperIter]+(nperIter*(l-1))) = CcrossValid{ii}(i,:);
        end
        
        nets{i,ii}([1:nperIter]+(nperIter*(l-1))) = Net{ii,i};
        
       end
        
   end
   
     end

	[minIter,minIterInd] = nanmin(data,[],3);

   for i = 1:nsub

	[~,BEST_PARAM] = nanmin(minIter(i,:));
	
	[~,BEST_NET_IND] = nanmin(squeeze(data(i,BEST_PARAM,:)));

	Fcv.bestnet{j}{i} = nets{i,BEST_PARAM}{BEST_NET_IND};
	Fcv.bestnetcorr{j}(i) = cdata(i,BEST_PARAM,BEST_NET_IND);


   end
   
   disp('Finding best networks')

   bestNetPerParam = cell(1,nsub);
  for i = 1:100
    ITER = 1;
    for ii = 1:100
        if ii ~= i
        	bestNetPerParam{i}{ITER} = nets{i,ii}{minIterInd(i,ii)};
        	ITER = ITER + 1;
        end
    end
  end	
	
%     display(['Getting model properties'])
% 	for i = 1:100
%         tic
% 		for ii = 1:100
%             NetProp = zeros(20,4);
%             NetProp_nets = nets{i,ii};
% 			for iii = 1:20
% 				A = zeros(100);
% 				b = NetProp_nets{iii};
% 				A(b) = 1;
% 				A = A + A';
% 				NetProp(iii,:) = NetworkProperties(A);
%             end
%             Fcv.NetProps{j}{i,ii} = NetProp;
%         end
%         toc
% 	end

% cdf_data = [];
% for i = 1:100
%     for ii = 1:100
%         NetProp_nets = nets{i,ii};
%         cdf_data_temp1 = cell(20,1);
%         cdf_data_temp2 = cell(20,1);
%         cdf_data_temp3 = cell(20,1);
%         cdf_data_temp4 = cell(20,1);
%         parfor iii = 1:20
% 				A = zeros(99);
% 				b = NetProp_nets{iii};
% 				A(b) = 1;
% 				A = A + A';   
%                 cdf_data_temp1{iii} = sum(A);
%                 cdf_data_temp2{iii} = betweenness_bin(A);
%                 cdf_data_temp3{iii} = clustering_coef_bu(A);
%                 cdf_data_temp4{iii} = Dists(triu(A,1) > 0);
%         end
%         cdf_data = [cdf_data; cdf_data_temp1 cdf_data_temp2 cdf_data_temp3 cdf_data_temp4];
%     end
% end
% 
% 
% CDFmeasures{j} = cdf_data;

        timesec = toc;
        disp(['Finished ',num2str(j),' in ',num2str(timesec),' seconds'])

   NETS{j} = nets;
   Fcv.SUBDATA{j} = data;
   Fcv.corr{j} = cdata;
   Fcv.minmin{j} = nanmin(nanmin(data,[],3),[],2)';
   Fcv.meanmean{j} = nanmean(nanmean(data,3),2)'; 
   Fcv.meanmin{j} = nanmean(nanmin(data,[],3),2)';
   Fcv.minmean{j} = nanmin(nanmean(data,3),[],2)'; 
   Fcv.MinParamNet{j} = bestNetPerParam;
   
   Fcv.corr_maxmax{j} = nanmax(nanmax(cdata,[],3),[],2)';
   Fcv.corr_meanmean{j} = nanmean(nanmean(cdata,3),2)'; 
   Fcv.corr_meanmax{j} = nanmean(nanmax(cdata,[],3),2)';
   Fcv.corr_maxmean{j} = nanmax(nanmean(cdata,3),[],2)'; 
end



% savename = ['CrossValid_',TYPE,'_indvdist_',num2str(INDVDIST),'_growth_',num2str(GROWTH)];
% save([savename,'.mat'],'Fcv')
% 
% if savenets
%    save([savename,'_NETS.mat'],'NETS','-v7.3') 
% end
