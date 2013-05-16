function result = compare(c1,c2, tolerance)

%tstoolbox/@core/compare
%   Syntax:
%     * compare(c1,c2, tolerance)
%
%   Input Arguments:
%     * c1,c2 core object of two signals
%     * tolerance tolerance of the signals's RMS value (default
%       tolerance=1e-6)
%
%   compare compare two signals whether they have equal values slight
%   differences due to rounding errors are ignored depending on the value
%   of tolerance when signals are found to be not equal, a zero is
%   returned.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

if nargin < 2	
	help(mfilename)
	return
end

if nargin < 3	
	tolerance = 1e-6;
end

if dlens(c1) ~= dlens(c2)
	result = 0;
else
	ref = data(rms(c1));
	diff = data(rms(c1-c2));
	if any((diff ./ ref) > tolerance)
		result = 0;
	else
		result = 1;
	end
end

