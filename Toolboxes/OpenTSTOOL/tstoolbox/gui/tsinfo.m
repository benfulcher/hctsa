function status = tsinfo(filename)

% status = tsinfo(filename) 
% display information about signal file 'filename'


trenner = '--------------------------------------------------------------------------------';

if exist(filename)
	sig = signal(filename);
	if isempty(sig)
		status = -1;
		warning(['Error opening file ' filename ]);
		return;
	else
		temp = dir(filename);
		filesize = num2str(temp.bytes);
		[path,Name,ext,ver] = fileparts(filename);
		yunt = unit(yaxis(sig));		
		string1 = { ...
		['Filename                        : ' filename] ; ...
		['Filesize                         : ' filesize ' bytes'] ; ...
		['Number of dimensions       : ' num2str(length(dlens(sig)))] ; ...
		['Length in each dimension   : ' num2str(dlens(sig))] ; ...
		['Total samples                  : ' num2str(prod(dlens(sig)))] ; ...
		trenner ; ...
 		['yname                           : ' name(yaxis(sig))] ; ...
		['yunit                              : '  name(yunt)] ; ...
		trenner };
		string2 = {};
% 		siz = size(sig.optparam);
% 		ende = siz(1);			
% 		for i=1:ende
% 			string2 = addcomment(string2, [sig.optparam{i}]);
% 		end
%		string2 = addcomment(string2, trenner);
% 		string3 = {};
% 	   for i=1:length(dlens(sig))
% 			xax = xaxis(sig, i);
% 		   string3 = addcomment(string3, ['Dimension : ' num2str(i)]);
% 		   string3 = addcomment(string3, ['Length : ' num2str(dlens(sig, i))]);
% 		   string3 = addcomment(string3, ['xfirst : ' num2str(first(xax))]);
% 		   string3 = addcomment(string3, ['xdelta : ' num2str(delta(xax))]);
% 		   %string3 = addcomment(string3, ['xdelta^(-1) : ' num2str(1/delta(xax))]);
% 		   string3 = addcomment(string3, ['xunit : ' label(xax)]);
% 		   string3 = addcomment(string3, ['xname : ' name(xax)]);
% 		   string3 = addcomment(string3, trenner);
% 	   end
% 
% 		string4 = { ['Comment :'] ; ''};		
% 		string = addcomment(string3, trenner);
% 		string5 = { trenner; ['History :'] ; ''};		
% 		string = cat(1, string1, string2, string3, string4, comment(sig), string5, history(sig));
% 		
		string = string1;
		handle = figure('Color', [1 1 1], ...
		'Resize','off', ...
		'NumberTitle','off', ...
		'KeyPressFcn','close(gcbf)', ...
		'MenuBar','none', ...
		'Name',['tsinfo on file ' Name]);
		pos = get(handle, 'Position');
		lbox = uicontrol('Parent',handle, ...
		'Units','points', ...
		'Position', [0.8 0.8 440 330], ...
		'BackgroundColor',[1 1 1], ...
		'String', string, ...
		'Style','listbox', ...
		'Tag','InfoListbox', ...
		'Value',1);
	end
else
	warndlg(['Can''t find file ' filename],'Info on File');
	status = -1;
	return;
end
