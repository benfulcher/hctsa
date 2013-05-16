function n = name(d)

%tstoolbox/@description/name
%   description/name Syntax:
%     * n = name(d)
%
%   Get signal's name
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

error(nargchk(1,1, nargin));

n = d.name;
