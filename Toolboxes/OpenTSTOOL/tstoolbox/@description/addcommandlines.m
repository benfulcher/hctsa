function d = commandlines(d, commandname, varargin)

%tstoolbox/@description/addcommandlines
%   adds new commandline to list of commands that have been applied to
%   that signal
%   example 1
%          addcommandlines(s, 's = spec2(s', 512, 'Hanning' )) will add 's
%          = spec2(s, 512, 'Hanning');' to the list of applied commands
%   example 2
%          len = 512; text = 'Hanning'; addcommandlines(s, 's = spec2(s',
%          512, 'Hanning' )) will add 's = spec2(s, 512, 'Hanning');' to
%          the list of applied commands
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt


commandline = commandname;

for i=1:nargin-2
	commandline = [commandline ', ' tostring(varargin{i})];
end

commandline = [commandline ');'];
d.commandlines = append(d.commandlines, commandline);

function out = tostring(in)
switch class(in)
	case 'char'
		out = ['''' in ''''];
	case 'double'
		out = num2str(in);
	otherwise
		error('Wrong type of argument(s) given');
end
