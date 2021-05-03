function lgd = VisualiseFitsLgd(group1,group2,group3,pos)

for j = 1:3
    
    if j == 1
        xtick_name = group1;
    elseif j == 2
        xtick_name = group2;
    elseif j == 3
        xtick_name = group3;
    end

if ~isempty(xtick_name)
    filesubind = strfind(xtick_name,'\n');    
    if ~isempty(filesubind)    
    LineEnds = [-2 filesubind-1 length(xtick_name)];
        for i = 1:length(LineEnds)-1
            tick_labels{i,j} = xtick_name(LineEnds(i)+3:LineEnds(i+1));            
        end
    else
       tick_labels{1,j} = xtick_name;  
    end
    
end

end

Ntypes2Plot = find([~isempty(group1) ~isempty(group2) ~isempty(group3)]);
tick_labels = tick_labels(:,Ntypes2Plot);

Nlines = size(tick_labels,1);

if Nlines == 1
    
    tickLabels = tick_labels;
    
else
    
    format_new_lines = ['%s',repmat('\\newline%s',1,Nlines-1),'\n'];
    
    %tickLabels = strtrim(sprintf(format_new_lines, tick_labels{:}));
    tickLabels = sprintf(format_new_lines, tick_labels{:});
end

lgd = axes('Position',pos);

    x1 = [-3:.1:3];
    x2 = x1 + 6.6;
    y1 = normpdf(x1,0,1);
    
    Nfakedata = 100;
    
    N2gen = ceil(y1*Nfakedata);
    FakeData1 = [];
    for i = 1:length(y1)
        FakeData1 = [FakeData1 ones(1,N2gen(i))*x1(i)];
    end

    FakeData2 = [];
    for i = 1:length(y1)
        FakeData2 = [FakeData2 ones(1,N2gen(i))*x2(i)];
    end
    
    for i = Ntypes2Plot
        %Violin(FakeData1,i,'ShowData',false,'ViolinColor',Colors2use(i,:),'EdgeColor',findAlphaColor(Colors2use(i,:),1),'ViolinAlpha',.5);
        V1 = Violin(FakeData2,i,'ShowData',false,'ViolinColor',[0 0 0],'EdgeColor',[0 0 0],'ViolinAlpha',0);
    if i == 1
        V1.EdgeColor = [0 0 0];
        V1.ViolinPlot.LineWidth = 2;
        V1.ViolinPlot.LineStyle = '-';
    elseif i == 2
        V1.EdgeColor = [0 0 0];
        V1.ViolinPlot.LineWidth = 2;
        V1.ViolinPlot.LineStyle = '--';
        
    elseif i == 3
        V1.EdgeColor = [0 0 0];
        V1.ViolinPlot.LineWidth = 2;
        V1.ViolinPlot.LineStyle = ':';
    end

    end
    

%     yticks([0 6.6])
%     ylim([-3.4 10])
%     xticks(1:5)
%     xlim([.5 5.5])
    box on

    xticks(Ntypes2Plot)
    xticklabels(tickLabels)
    yticks([])
    
    %xticklabels({'Spatial','Homophily','Clustering','Degree','Communicability'})
    set(gca,'FontSize',14)
    
end