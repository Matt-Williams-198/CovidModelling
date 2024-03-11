
for i = 1:size(Splot,1)
    plot(x,squeeze(StoredDataI(1,i,:)),'b');
    hold on
    plot(x,squeeze(StoredDataI(10,i,:)),'r');
    plot(x,squeeze(StoredDataI(20,i,:)),'k');
    xlabel('Spatial Dimension','FontSize',18,'interpreter','latex')
    ylabel('$S$, $I$, and $R$','FontSize',18,'interpreter','latex')
    set(gca,'fontsize',17)
    ylim = max(Splot);
    drawnow; 
    pause(0.05)
    hold off
end
