x = linspace(0,ModelLength,SpatialDiscretization);
for i = 1:length(SIRS)
    plot(x,SIRS(i,:),'b');
    hold on
    plot(x,SIRI(i,:),'r');
    plot(x,SIRR(i,:),'k');
    xlabel('Spatial Dimension','FontSize',18,'interpreter','latex')
    ylabel('$S$, $I$, and $R$','FontSize',18,'interpreter','latex')
    set(gca,'fontsize',17)
    drawnow; 
    pause(1)
    hold off
end