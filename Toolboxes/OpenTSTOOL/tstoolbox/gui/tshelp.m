function tshelp(funcname)

% tshelp(functionname)
%
% displays the help for a given function <functionname>
% in a separate window
%
% cmerk 1998

hlptext = help(['signal/' funcname]);

if ~isempty(hlptext)
	title = ['Help for ' funcname];
	lineNo = 12;
	prompt = {''};
	answer  = inputdialog(prompt,title,lineNo, {hlptext});
else
	warndlg(['No help for ' funcname 'available']);
end
