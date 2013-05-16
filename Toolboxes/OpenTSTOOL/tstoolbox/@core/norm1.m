function cout = norm1(cin, low, upp)

%tstoolbox/@core/norm1
%   Syntax:
%     * cout = norm1(cin,low,upp)
%
%   Input Arguments:
%     * cin - core object
%     * low - column number
%     * upp - column number
%
%   normalize each single column of a the core object to be within
%   [low,upp]
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

N = dlens(cin,1);	
points = data(cin);
points = points - repmat(min(points), N,1);
points = points ./ repmat(max(points)/(upp-low), N,1);
if low ~= 0
	points = points + low;
end
cout = core(points);
