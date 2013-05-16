function v = spacing(a, len)

%tstoolbox/@achse/spacing
%   Syntax:
%     * v = values(a, len)
%
%   Returns spacing values for linear, logarithmic or arbitary spacing in
%   case of lin. or log. spacing. len values are returned. In case of
%   arbitary spacing, all stored values are returned.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt


if nargin < 1, help(mfilename), return, end
	
if isa(a, 'achse')
	switch a.resolution
		case 'linear'
			if nargin < 2, help(mfilename), return, end
			v = a.first:a.delta:(a.first+(len-1)*a.delta);
		case 'logarithmic'
			if nargin < 2, help(mfilename), return, end
			v = a.first * ( a.delta .^ (0:len-1) );
		case 'arbitrary'
			v = a.values;
	end
else
	error('need achse object');
end


