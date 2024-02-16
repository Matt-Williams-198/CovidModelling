figure
h = line(nan, nan, 'Color', 'b');
% Loop over the data to animate the line
for k = 1:length(squeeze(SystemIterationsS(1,1,:)))
    set(h, 'XData', SpatialVector, 'YData', squeeze(SystemIterationsS(1,k,:)));
    drawnow;
end
