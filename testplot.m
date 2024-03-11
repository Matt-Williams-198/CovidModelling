figure
SLine = line(nan, nan, 'Color', 'b');
ILine = line(nan, nan, 'Color', 'g');
RLine = line(nan, nan, 'Color', 'r');
set(SLine, 'XData', SpatialVector(50:950), 'YData', Sarray(k,50:950));
set(ILine, 'XData', SpatialVector(50:950), 'YData', Iarray(k,50:950));
set(RLine, 'XData', SpatialVector(50:950), 'YData', Rarray(k,50:950));
drawnow;
pause(1);
for k = 1:length(Sarray(:,1))
    set(SLine, 'XData', SpatialVector(50:950), 'YData', Sarray(k,50:950));
    set(ILine, 'XData', SpatialVector(50:950), 'YData', Iarray(k,50:950));
    set(RLine, 'XData', SpatialVector(50:950), 'YData', Rarray(k,50:950));
    drawnow;
    pause(1)
end
