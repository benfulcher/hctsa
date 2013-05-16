function a = getaxis(s, dim)

%tstoolbox/@signal/getaxis
%   Syntax:
%     * a = getaxis(s, dim)
%
%   Get one of the currend xaxes.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt
error(nargchk(1,2, nargin));

if nargin < 1
	help(mfilename)
	return
end

if nargin < 2
	a = s.xaxes{1};
else
	a = s.xaxes{dim};	% FIXME - no range checking
end
