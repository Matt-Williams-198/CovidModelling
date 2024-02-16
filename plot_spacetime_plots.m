function plot_spacetime_plots(name, t_max, st_arrayS, st_arrayI)
    % Plot the space-time plots
    % Plot S 

    figure('Position', [10, 10, 340, 250]);
    imagesc([0, L], [0, t_max], st_arrayS);
    colormap('hot');
    caxis([0.001, 1.0]);
    colorbar;
    xlabel('$x$', 'Interpreter', 'latex');
    ylabel('$t$', 'Interpreter', 'latex');
    print('-dpdf', strcat(name, '/S.pdf'));

    % Plot I 
    figure('Position', [10, 10, 340, 250]);
    imagesc([0, L], [0, t_max], st_arrayI);
    colormap('hot');
    caxis([0.001, 1.0]);
    colorbar;
    xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 10);
    ylabel('$t$', 'Interpreter', 'latex', 'FontSize', 10);
    print('-dpdf', strcat(name, '/I.pdf'));
end