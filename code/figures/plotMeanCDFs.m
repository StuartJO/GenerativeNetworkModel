function plotMeanCDFs(SubNets,NETS,Dists)

N = length(SubNets);

for j = 1:length(NETS)
    
    if size(NETS{j},1) ~= size(NETS{j},2)
        A = zeros(N);
        b = NETS{j};
        A(b) = 1;
        NETS{j} = A + A';
    end
    
    SubNet = double(SubNets{j}>0);
    
    subdata{j,1} = sum(SubNet);
    subdata{j,3} = betweenness_bin(SubNet);
    subdata{j,2} = clustering_coef_bu(SubNet);
    subdata{j,4} = Dists(triu(SubNet,1) > 0);
    
    data{j,1} = sum(NETS{j});
    data{j,3} = betweenness_bin(NETS{j});
    data{j,2} = clustering_coef_bu(NETS{j});
    data{j,4} = Dists(triu(NETS{j},1) > 0);
    
end

Colors = lines(4);

figure('Position',[506 536 1539 419])

x_labels = {'Degree, \itk','Clustering, \itc','Betweenness, \itb','Edge Length, \ite'};
y_labels = {'CDF(\itk)','CDF(\itb)','CDF(\itc)','CDF(\ite)'};

for i = 1:4

ax.plot(i) = subplot(1,4,i);


    [f,x] = mean_ecdf(data(:,i));
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

    [f,x] = mean_ecdf(subdata(:,i));
    plot(x,f,'Color',[0 0 0],'LineWidth',2)
end

ax.plot(1).Position([2 4]) = [0.2 0.6283];
ax.plot(2).Position([2 4]) = [0.2 0.6283];
ax.plot(3).Position([2 4]) = [0.2 0.6283];
ax.plot(4).Position([2 4]) = [0.2 0.6283];
% l = legend([NAMES(2:end) NAMES(1)],'Location','southeastoutside','Orientation','horizontal','NumColumns',6);
% 
% l.Position(1) = 0.5 - l.Position(3)/2; 
% l.Position(2) = 0.02;

% annotation('textbox',...
%     [0.13 0.863961813842482 0.775178687459389 0.124105011933174],...
%     'String',TITLE,...
%     'FontSize',24,...
%     'FitBoxToText','off','EdgeColor','none','HorizontalAlignment','center');