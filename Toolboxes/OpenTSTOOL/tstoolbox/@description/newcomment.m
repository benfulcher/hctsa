function d = newcomment(d, c)

%tstoolbox/@description/newcomment
%   Syntax:
%     * d = newcomment(d, string)
%     * d = newcomment(d, list)
%
%   Replace old comment with new comment
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(2,2);	

d.comment = list(c);	% c may have every type that is accepted by list()

