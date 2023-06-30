function [] = statAnalySSM(SSM)
% show statistical results of SSM


figure
scatter3(SSM.bVecs(:,1),SSM.bVecs(:,2),SSM.bVecs(:,3));
xlabel('component 1')
ylabel('component 2')
zlabel('component 3')
title('position of datapoints on the first 3 PCs')

prop=zeros(1,length(SSM.exp));
for i=1:length(SSM.exp)
    prop(i)=sum(SSM.exp(1:i));
end
figure
% pareto(SSM.eVals);
yyaxis left
bar(SSM.eVals)
ylabel('eigenvalues')
yyaxis right
plot(prop)
ylabel('percentage')
xlabel('componnents')
xticks(1:1:length(SSM.exp));
xt = xticks(gca);
xt2 = compose('PC%d',xt);
set(gca,'xtickLabel',xt2);
ylim([0 100])
title('eigenvalue and percentage distribution of PCs')
end