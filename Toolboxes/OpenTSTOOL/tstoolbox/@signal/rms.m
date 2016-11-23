function rs = rms(s)

%tstoolbox/@signal/rms
%   Syntax:
%     * rs = rms(s)
%
%   Calculate root mean square value for signal along dimension 1.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(1,1);

c = rms(s.core); 	% call real working routine for parent core object
rs = signal(c, s);

if ndim(rs) > 0
	rs.xaxes = s.xaxes(2:end);
end

rs = addhistory(rs,  ['Calculated root mean square value along dimension 1'] );
rs = addcommandlines(rs, 's = rms(s');
