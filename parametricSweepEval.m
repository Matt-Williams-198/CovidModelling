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
t = linspace(0,tfinal,100);
plot(t,squeeze(PDIofAlpha(1,:)),'b');
hold on
plot(t,squeeze(PDIofAlpha(10,:)),'r');
plot(t,squeeze(PDIofAlpha(20,:)),'k');
xlabel('Time','FontSize',18,'interpreter','latex')
ylabel('Mean P.D. of $I$','FontSize',18,'interpreter','latex')
legend('alpha = 0.005','alpha = 0.05','alpha = 0.1')
set(gca,'fontsize',17)
ylim = max(Splot);
hold off
