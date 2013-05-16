function linefit
zoom off
k = waitforbuttonpress;
point1 = get(gca,'CurrentPoint');    % button down detected
finalRect = rbbox;                   % return Figure units
point2 = get(gca,'CurrentPoint');    % button up detected
point1 = point1(1,1:2);              % extract x and y
point2 = point2(1,1:2);
p1 = min(point1,point2);             % calculate locations
offset = abs(point1-point2);         % and dimensions
x1 = p1(1);
x2 = p1(1)+offset(1);

zoom on
linehandle = findobj(gca, 'Type', 'line');
if ~isempty(linehandle)
    xdata = get(linehandle(1), 'XData');
    ydata = get(linehandle(1), 'YData');
    ind = find(xdata >=x1  & xdata <= x2);     
    [p,s] = polyfit(xdata(ind), ydata(ind), 1);
    if p(2) >= 0
        msgbox(['y = ' num2str(p(1)) ' * x + ' num2str(p(2))], 'Line fit');
    else
        msgbox(['y = ' num2str(p(1)) ' * x - ' num2str(abs(p(2)))], 'Line fit');
    end
end
