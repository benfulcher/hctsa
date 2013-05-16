function cout = norm2(cin)

%tstoolbox/@core/norm2
%   Syntax:
%     * cout = norm2(cin)
%
%   Input Arguments:
%     * cin - core object
%
%   normalize signal by removing it's mean and dividing by the standard
%   deviation
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

N = dlens(cin,1);	
points = data(cin);
points = points - repmat(mean(points), N,1);
points = points ./ repmat(std(points),N,1);
cout = core(points);
