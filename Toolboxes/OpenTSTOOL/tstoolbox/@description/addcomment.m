function d = addcomment(d, text)

%tstoolbox/@description/addcomment
%   adds new comment to current list of comments
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

error(nargchk(1,2,nargin));

d.comment = append(d.comment, text);
