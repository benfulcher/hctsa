function n = name(d)

%tstoolbox/@description/name
%   description/name Syntax:
%     * n = name(d)
%
%   Get signal's name
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(1,1);

n = d.name;
