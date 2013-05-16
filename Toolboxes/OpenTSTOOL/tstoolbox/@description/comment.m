function c = comment(d)

%tstoolbox/@description/comment
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

error(nargchk(1,1,nargin));

c = cellstr(d.comment);

