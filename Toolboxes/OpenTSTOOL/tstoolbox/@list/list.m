function l = list(varargin)

%tstoolbox/@list/list
%   Syntax:
%     * l = list
%       creates empty list
%     * l = list('Hello world')
%       create list with one entry, 'Hello world'
%     * l = list('Hello', 'My' , 'World')
%       create list with three entries
%     * l = list('Hello', 'My' , 'World')
%       create list with three entries
%
%   An object of type list contains a list of strings.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt


if nargin == 0 			% create an empty list
	l.data = {};		% implementation as one-dimensional cell array (columnwise)
	l.len = 0;
	l = class(l, 'list');
else					% create list from arguments
	if iscellstr(varargin)
		l = list(varargin);
	elseif ischar(varargin{1})
	    l = list(cellstr(varargin{1}));
	elseif isa(varargin{1}, 'list')
		l = varargin{1};
	elseif iscellstr(varargin{1})
		l.data = varargin{1}(:);
		l.len = length(varargin{1});
		l = class(l, 'list');
	else
		error('Wrong type of argument(s) given');
	end
end

