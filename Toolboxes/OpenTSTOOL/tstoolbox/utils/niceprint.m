function niceprint(fileName, varargin)

% niceprint(fileName, fontSize, format,  figHandle)
%
% Prepare a figure (default is current figure,
% otherwise with handle <figHandle>) for printing to a file
% named <fileName> with fontsize <fontSize> as color eps level II
% (default), using some default settings for a nice appearance
%
% Format may be :
%
% eps - color encapsulated postscript level II
% ps - color postscript level II
% tiff - tiff image file
%
% and for use with the zbuffer mode :
%
% zeps - color encapsulated postscript level II
% zps - color postscript level II
% ztiff - tiff image file

narginchk(1,4);

if nargin < 4
	figHandle = gcf;
else
	figHandle = varargin{3};
end
if nargin < 3
	format = 'eps';
else
	format = varargin{2};
end
if nargin > 1
	fontSize = varargin{1};
end

figure(figHandle)
axHandle = gca;

% Enable usage of Tex characters in labels
set(figHandle, 'DefaultTextInterpreter', 'Tex')

% Set axis properties
set(axHandle, 'PlotBoxAspectRatioMode', 'auto')
% set(axHandle, 'PlotBoxAspectRatio', [px py pz])
set(axHandle, 'Units', 'normalized');
Width = 0.70;
XYRatio = 0.7;	% Height/Width
Height = XYRatio * Width;
set(axHandle, 'Position', [(1-Width)/2 (1-Height)/2 Width Height])

xl = get(gca, 'xlabel');
set(xl, 'Units', 'normalized')
set(xl, 'Interpreter', 'Tex');
position = get(xl, 'Position');
set(xl, 'Position', [0.4954   -0.18 0]);

yl = get(gca, 'ylabel');
set(yl, 'Units', 'normalized')
set(yl, 'Interpreter', 'Tex');
position = get(yl, 'Position');
set(yl, 'Position', [ -0.11 0.4942 0]);

tl = get(gca, 'title');
set(tl, 'Units', 'normalized')
set(tl, 'Interpreter', 'Tex');
position = get(tl, 'Position');
set(tl, 'Position', [0.4954 1.00 0]);

% Chnage FontSize of all axes and text objects
set(0, 'ShowHiddenHandles', 'on')

set(axHandle, 'FontUnits', 'points');
axesFontSize = get(axHandle, 'FontSize');
scaleFactor = fontSize / axesFontSize;

handles = [findall(gcf, 'Type', 'text'); findall(gcf, 'Type', 'axes')];

for h=handles(:)'
	set(h, 'FontSize', get(h, 'FontSize')*scaleFactor);
end

%
% for c=childs(:)'		% '
% 	if strcmp(get(c, 'Type'), 'text')  | strcmp(get(c, 'Type'), 'axes')
% 		set(c, 'FontSize', fontSize)
% 	end
% % 	if strcmp(get(c, 'Type'), 'line')
% % 		set(c, 'MarkerSize', 36)
% % 		set(c, 'Marker', '.')
% % 	end
% end
%

set(0, 'ShowHiddenHandles', 'off')

% Set global printing properties
set(figHandle, 'Color', [1 1 1]);
%set(figHandle, 'Renderer', 'zbuffer');
set(figHandle, 'PaperType', 'A4');
%set(figHandle, 'PaperOrientation', 'landscape');
set(figHandle, 'PaperOrientation', 'portrait');
set(figHandle, 'PaperUnits', 'normalized');

% Set Paper position
%set(figHandle, 'PaperPositionMode', 'auto');
% or
set(figHandle, 'PaperPositionMode', 'manual');
Width = 0.7;
XYRatio = 1;	% Height/Width
Height = XYRatio * Width;
set(figHandle, 'PaperPosition', [(1-Width)/2 (1-Height)/2 Width Height]);

resolution = 80;
switch format
	case 'ps'
		print(['-f' num2str(figHandle)], '-noui', '-dpsc2', '-painters', fileName)
	case 'tiff'
		print(['-f' num2str(figHandle)], '-noui', '-dtiff', '-painters', fileName)
 	case 'zeps'
		print(['-f' num2str(figHandle)], ['-r' num2str(resolution)], '-noui', '-depsc2', '-zbuffer', fileName)
	case 'zps'
		print(['-f' num2str(figHandle)], ['-r' num2str(resolution)], '-noui', '-dpsc2', '-zbuffer', fileName)
	case 'ztiff'
		print(['-f' num2str(figHandle)], ['-r' num2str(resolution)] , '-noui', '-dtiff', '-zbuffer', fileName)
	otherwise
		print(['-f' num2str(figHandle)], '-noui', '-depsc2', '-painters', fileName)
end

%unix(['gv ' fileName ' &']);
unix(['xpsview ' fileName ' &']);
