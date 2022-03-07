function plotCDFs(subdata,MdlData,Type,Colors)

% for j = 1:length(SubNets)
%     
%     SubNet = double(SubNets{j}>0);
%     
%     subdata{j,1} = sum(SubNet);
%     subdata{j,3} = betweenness_bin(SubNet);
%     subdata{j,2} = clustering_coef_bu(SubNet);
%     subdata{j,4} = Dists(triu(SubNet,1) > 0);    
%     
% end

% figure('Position',[506 536 1539 419])

x_labels = {'Degree, \itk','Clustering, \itc','Betweenness, \itb','Node mean distance, \ite'};
y_labels = {'CDF(\itk)','CDF(\itb)','CDF(\itc)','CDF(\ite)'};

for i = Type

% ax.plot(i) = subplot(1,4,i);

    if size(MdlData{i},1) == 1
        [f,x] = ecdf(MdlData{i});    
    else
        [f,x] = mean_ecdf(MdlData{i});
    end
    plot(x,f,'Color',Colors(i,:),'LineWidth',2)
hold on
box on
% for j = 1:length(NETS)
% [f,x] = ecdf(data{j,i});
% 
% plot(x,f,'Color',make_alpha_rgb(Colors(i,:),.5),'LineWidth',.5)
% end


    xlabel(x_labels{i})
    ylabel(y_labels{i})
    
    %legend('Empirical','Additive','Multiplicative','Location','southeast')

    set(gca,'FontSize',16)
    if size(MdlData{i},1) == 1
        [f,x] = ecdf(subdata{i});    
    else
        [f,x] = mean_ecdf(subdata{i});
    end
    plot(x,f,'Color',[0 0 0],'LineWidth',2)
end