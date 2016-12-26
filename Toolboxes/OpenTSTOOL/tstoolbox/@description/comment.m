function c = comment(d)

%tstoolbox/@description/comment
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(1,1);

c = cellstr(d.comment);

