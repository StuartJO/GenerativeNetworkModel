function CompareEdgeLengtheCDF(A_dist,MultNets,Add2Nets,Add3Nets,adjs)

% for i = 1:100; Net = Fcv_mult2_static.bestnet{3}{i}; D = Net.*dists{19}; Du = triu(D,1); DATA{i} = Du(Du~=0); end
% for i = 1:100; Net = Fcv_add2_static.bestnet{3}{i}; A = zeros(100); A(Net) = 1; A = A + A'; Net = A; D = Net.*dists{19}; Du = triu(D,1); DATA2{i} = Du(Du~=0); end
% for i = 1:100; Net = Fcv_add3_static.bestnet{3}{i}; A = zeros(100); A(Net) = 1; A = A + A'; Net = A; D = Net.*dists{19}; Du = triu(D,1); DATA4{i} = Du(Du~=0); end
% for i = 1:100; Net = double(adjs{i}>0); D = Net.*dists{19}; Du = triu(D,1); DATA3{i} = Du(Du~=0); end

for i = 1:length(MultNets)
    Du = triu(MultNets{i}.*A_dist,1);
    MultData{i} = Du(Du~=0);
end
for i = 1:length(Add2Nets)
    Du = triu(Add2Nets{i}.*A_dist,1);
    Add2Data{i} = Du(Du~=0);
end
for i = 1:length(Add3Nets)
    Du = triu(Add3Nets{i}.*A_dist,1);
    Add3Data{i} = Du(Du~=0);
end
for i = 1:length(adjs)
    Du = triu(double(adjs{i}>0).*A_dist,1);
    AdjData{i} = Du(Du~=0);
end

for i = 1:length(MultData); [f(i,:),xi] = ksdensity(MultData{i},0:160);  end
for i = 1:length(Add2Data); [f1(i,:),xi] = ksdensity(Add2Data{i},0:160);  end
for i = 1:length(Add3Data); [f2(i,:),xi] = ksdensity(Add3Data{i},0:160); end
for i = 1:length(AdjData); [f3(i,:),xi] = ksdensity(AdjData{i},0:160); end

plot(xi,mean(f),'LineWidth',2)
hold on
plot(xi,mean(f1),'LineWidth',2)
plot(xi,mean(f2),'LineWidth',2)

plot(xi,mean(f3),'Color',[0 0 0],'LineWidth',2)
ylabel('Proportion of connections')
xlabel('Connection length')
legend({'Multiplicative','Additive (no \gamma)','Additive','Empirical'})
set(gca,'FontSize',24)