function rs = center(s)

%tstoolbox/@signal/center
%   Syntax:
%     * center(s)
%
%   Center signal by removing it's mean.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(1,1);
if nargin  < 1, help(mfilename); end

N = dlens(s,1);
points = data(s);
c = core(points - repmat(mean(points), N,1));
rs = signal(c,s);

rs = addhistory(rs, ['Centered signal around zero']);
rs = addcommandlines(rs, 's = center(s');
