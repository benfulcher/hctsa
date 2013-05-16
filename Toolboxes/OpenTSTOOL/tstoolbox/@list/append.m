function l = append(l, argument)

%tstoolbox/@list/append
%   Syntax:
%     * list = append(list, string)
%     * list = append(list, list)
%
%   Add string(s) to existing list.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

l.data = l.data(:); 		% ensure compatibility with previous implementation

if ischar(argument)
    l.data = [l.data; {argument}];
	l.len = l.len + 1;
elseif isa(argument, 'list')
	l.len = l.len + argument.len;
	l.data = [l.data; argument.data];
elseif iscellstr(argument)
	l = append(l, list(argument));
else
	error('Wrong type of argument(s) given');
end
