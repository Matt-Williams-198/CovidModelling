PDSofAlpha = zeros(size(StoredDataS,1),size(StoredDataS,2));
PDIofAlpha = zeros(size(StoredDataS,1),size(StoredDataS,2));
PDRofAlpha = zeros(size(StoredDataS,1),size(StoredDataS,2));
for i = 1:size(StoredDataS,1)
    for j = 1:size(StoredDataS,2)
        PDSofAlpha(i,j) = sum(StoredDataS(i,j,:))/ModelLength;
        PDIofAlpha(i,j) = sum(StoredDataI(i,j,:))/ModelLength;
        PDRofAlpha(i,j) = sum(StoredDataR(i,j,:))/ModelLength;
    end
end
alphaOutputs = zeros(size(PDIofAlpha,1),3);
for i = 1 : size(PDIofAlpha)
    alphaOutputs(i,1) = sum(PDSofAlpha(:,i));
    alphaOutputs(i,2) = sum(PDIofAlpha(:,i));
    alphaOutputs(i,3) = sum(PDRofAlpha(:,i));
end
t = linspace(0,tfinal,100);
plot(t,squeeze(PDIofAlpha(1,:)),'b');
hold on
plot(t,squeeze(PDIofAlpha(10,:)),'r');
plot(t,squeeze(PDIofAlpha(20,:)),'k');
title('Population Density over time against Alpha')
xlabel('Time','FontSize',18,'interpreter','latex')
ylabel('Mean P.D. of $I$','FontSize',18,'interpreter','latex')
legend('alpha = 0.005','alpha = 0.05','alpha = 0.1')
set(gca,'fontsize',17)
ylim = max(Splot);
hold off

TotalPop = sum(alphaOutputs(1,:));
alphaOutputs = alphaOutputs *10./TotalPop;
localx = linspace(0.1,1,20);

figure
plot(localx,alphaOutputs(:,2));
hold on
title('Total population against Alpha')
xlabel('Alpha')
ylabel('Total Infected Population')
plot(localx,alphaOutputs(:,1));
plot(localx,alphaOutputs(:,3));
legend('Infected', 'Susceptible','Recovered')
