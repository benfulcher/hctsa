function tsprint(figHandle, fileName, fontSize)

% tsprint(fileName, fontSize)
% tsprint(figHandle, fileName, fontSize)
%
% Prepare a figure with handle <figHandle> for printing to a file
% named <fileName> with fontsize <fontSize> as color postscript level II,
% using some default settings for a nice appearance

if nargin < 3
	figHandle = gcf;
end


fileName = fullfile(pwd, fileName);

set(0, 'ShowHiddenHandles', 'on')
childs = findobj(figHandle);

set(figHandle, 'DefaultTextInterpreter', 'Tex')
set(gca, 'PlotBoxAspectRatioMode', 'auto')
set(gca, 'Units', 'normalized');
set(gca, 'Position', [0.15 0.15 0.70 0.70])

xl = get(gca, 'xlabel');
set(xl, 'Units', 'normalized')
set(xl, 'Interpreter', 'Tex');
position = get(xl, 'Position');
set(xl, 'Position', [0.4954   -0.14 0]);

yl = get(gca, 'ylabel');
set(yl, 'Units', 'normalized')
set(yl, 'Interpreter', 'Tex');
position = get(yl, 'Position');
set(yl, 'Position', [ -0.11 0.4942 0]);

tl = get(gca, 'title');
set(tl, 'Units', 'normalized')
set(tl, 'Interpreter', 'Tex');
position = get(tl, 'Position');
set(tl, 'Position', [0.4954 1.0646 0]);

% FontSize auf einheitlichen Wert setzen

for c=childs(:)'		% '
	if strcmp(get(c, 'Type'), 'text')  | strcmp(get(c, 'Type'), 'axes')
		set(c, 'FontSize', fontSize)
	end
% 	if strcmp(get(c, 'Type'), 'line')
% 		set(c, 'MarkerSize', 36)
% 		set(c, 'Marker', '.')
% 	end
end

%set(figHandle, 'Renderer', 'zbuffer')
set(figHandle, 'PaperType', 'A4')
set(figHandle, 'PaperPositionMode', 'manual')
%set(figHandle, 'PaperPositionMode', 'auto')
set(figHandle, 'PaperOrientation', 'landscape')
%set(figHandle, 'PaperOrientation', 'portrait')
set(figHandle, 'PaperUnits', 'normalized')
set(figHandle, 'PaperPosition', [0.14 0.1 0.72 0.5])

print(['-f' num2str(figHandle)], '-noui', '-depsc2', [fileName])
%print(['-f' num2str(figHandle)], '-r50', '-noui', '-dtiff', '-zbuffer', [fileName '.tiff'])
%print(['-f' num2str(figHandle)], '-r100', '-noui', '-depsc', '-zbuffer', fileName)
%print(['-f' num2str(figHandle)], '-r100', '-noui', '-dtiff', '-zbuffer', fileName)
unix(['xpsview ' fileName ' &']);
set(0, 'ShowHiddenHandles', 'off')

