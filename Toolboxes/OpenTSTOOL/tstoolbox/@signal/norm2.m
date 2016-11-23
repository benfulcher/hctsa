function rs = norm2(s)

%tstoolbox/@signal/norm2
%   Syntax:
%     * rs=norm2(s)
%
%   Normalize signal by removing it's mean and dividing by the standard
%   deviation.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(1,1);

c = norm2(s.core);
rs = signal(c,s);

rs = addhistory(rs, ['(norm2) Centered signal and divided it by it''s standard deviation']);
rs = addcommandlines(rs, 's = norm2(s'); 

