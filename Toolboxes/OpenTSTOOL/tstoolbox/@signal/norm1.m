function rs = norm1(s, low, upp)

%tstoolbox/@signal/norm1
%   Syntax:
%     * rs=norm1(s) => low=0 , upp=1
%     * rs=norm1(s, low) => upp=1
%     * rs=norm1(s, low, upp)
%
%   Scale and move signal values to be within [low,upp].
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(1,3);

if nargin < 2
	low = 0;
end
if nargin < 3
	upp = 1;
end

c = norm1(s.core, low, upp);
rs = signal(c,s);

rs = addhistory(rs, ['(norm1) Transformed signal to be within [' ...
     num2str(low) ',' num2str(upp) ']']);
rs = addcommandlines(rs, 's = norm1(s');

