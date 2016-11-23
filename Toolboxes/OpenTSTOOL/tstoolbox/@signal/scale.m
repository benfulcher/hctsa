function rs = scale(s, f)

%tstoolbox/@signal/scale
%   Syntax:
%     * scale(signal, factor)
%
%   Scale signal by factor f.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(2,2);

c = core(data(s)*f);
rs = signal(c,s);

rs = addhistory(rs,  ['Scaled signal by factor ' num2str(f)]);
rs = addcommandlines(rs, 's = scale(s', f);
