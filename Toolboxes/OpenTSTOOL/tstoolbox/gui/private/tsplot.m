function tsplot(filename, varargin)


if nargin < 2
	mode = 'large'; 	% im Modus 'large' wird eine eigene Figure gestartet
else
	mode = varargin{1}; % im Modus 'small' wird in das Preview-Areal des tstool geplottet
end

if isunix		% use greater fonts on Unix workstations
	if strcmp(mode, 'small')
		fontsize = 14;
	else
		fontsize = 16;
	end
else
	if strcmp(mode, 'small')
		fontsize = 6;
	else
		fontsize = 9;
	end
end
if exist(filename, 'file')
	sig = signal(filename);
	dlen = dlens(sig);
	if  (ndim(sig)==1) & (dlen(1) < 5)
		[path,name,ext,ver] = fileparts(filename);
		titel = ['Values of ' name];
		lineNo = dlen(1);
		prompt  = 'Value(s) :';
		defaults = { num2str(data(sig)) };
		answer  = inputdlg(prompt,titel,lineNo,defaults);	% discard answer
	else
		if strcmp(mode, 'small')
			fhandle = findobj('Tag', 'TSTOOL' );
			if strcmp(plothint(sig),'subplotgraph') & (dlens(sig,2)>1)
			  sig=cut(sig,2,1,1);
			end
		else
			fhandle = newplotwin(filename);			
        end			

	view(sig, fontsize, fhandle);		
		if strcmp(mode, 'small')
			zoom off
			rotate3d off
		else
			set(fhandle, 'Visible', 'on');
		end
	end
else
	error('File not found');
end 								
	
		


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fhandle = newplotwin(filename)

% Create new figure window with own parameters
% Filename is only needed for the title of the window
scrsz = get(0,'ScreenSize');   % get Screensize

[path,name,ext,ver] = fileparts(filename);
titlestring = [name ext];
indices = find(titlestring=='_');
titlestring(find(titlestring=='_')) = ' ';  % remove all underscores

fhandle = figure(... % 'Position',[400 400 scrsz(3)/3.5 scrsz(4)/3] , ...
	'NumberTitle','off', ...
	'Name',titlestring, ...
	'Visible', 'off', ...
	'Units', 'normalized', ...
	'Color', [1 1 1]);

% b = uicontrol('Parent', fhandle, ...
% 	'Units','normalized', ...
% 	'Position',[0.050 0.01 0.92 0.035], ...
% 	'Callback', 'pos=get(gcbo,''Value'');limits=zoom(''getlimits'');axen=get(gcbf,''currentaxes'');XLim=get(axen,''XLim'');dx=XLim(2)-XLim(1);fl=limits(2)-limits(1)-dx;set(axen,''XLim'',[limits(1)+pos*fl limits(1)+pos*fl+dx])', ...
% 	'Style','slider');
% 	
% % b = axes('Parent', fhandle, ...
% % 	'Box','on', ...
% % 	'CameraUpVector',[0 1 0], ...
% % 	'Color',[1 1 1], ...
% % 	'Position',[0.050 0.20 0.92 0.70]);	
% 	

b =  uimenu('Parent',fhandle, ...
	'Label','Zoom');
	
c = uimenu('Parent',b, ...
	'Callback', 'limits=zoom(''getlimits'');axen=get(gcbf,''currentaxes'');YLim=get(axen,''YLim'');dy=YLim(2)-YLim(1);set(axen,''YLim'',[max(YLim(1)-(dy/2),limits(3)) min(YLim(2)+(dy/2),limits(4))])', ...
	'Label','Y out');
c = uimenu('Parent',b, ...
	'Callback', 'limits=zoom(''getlimits'');axen=get(gcbf,''currentaxes'');XLim=get(axen,''XLim'');dx=XLim(2)-XLim(1);set(axen,''XLim'',[max(XLim(1)-(dx/2),limits(1)) min(XLim(2)+(dx/2),limits(2))])', ...
	'Label','X out');
	
b = uimenu('Parent',fhandle, ...
	'Label','Fit');
	c = uimenu('Parent',b, ...
	'Callback','linefit', ...
	'Label','Line');

function p = get_currentpoint(ax)

p = get(ax,'currentpoint'); p = p(1,1:2);
if strcmp(get(ax,'XScale'),'log'),
  p(1) = log10(p(1));
end
if strcmp(get(ax,'YScale'),'log'),
  p(2) = log10(p(2));
end


function set_pos(ax)

pos_text=findobj('Tag','PositionsText');
if exist('pos_text')
	disp(pos_text);
	delete(pos_text);
else     
	disp('Kein pos_text');
end
position=get_currentpoint(ax)
pos_text=text(position(1),position(2),['<= (x=',num2str(position(1)),' y=',num2str(position(2)),' )'],'Tag','PositionsText');




