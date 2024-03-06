
for i = 1:size(Splot,1)
    plot(x,Splot(:,i),'b');
    hold on
    plot(x,Iplot(:,i),'r');
    plot(x,Rplot(:,i),'k');
    xlabel('Spatial Dimension','FontSize',18,'interpreter','latex')
    ylabel('$S$, $I$, and $R$','FontSize',18,'interpreter','latex')
    set(gca,'fontsize',17)
    ylim = max(Splot);
    drawnow; 
    pause(1)
    hold off
end