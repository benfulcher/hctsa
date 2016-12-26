function d = merge(d1, d2)

%tstoolbox/@description/merge
%   Syntax:
%     * d = merge(d1, d2)
%
%   merge two descriptions
%   Most items are taken from first description. History is taken from
%   both descriptions. This function may be useful when writing binary
%   operators for class signal
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(2,2);

d = d1; 	% copy all fields from first description;
d.history = list('{');
d.history = TSTL_append(d.history, d1.history);
d.history = TSTL_append(d.history, '}');
d.history = TSTL_append(d.history, '{');
d.history = TSTL_append(d.history, d2.history);
d.history = TSTL_append(d.history, '}');







