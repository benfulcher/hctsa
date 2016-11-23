function [D2, s] = takens_estimator(s, n, range, past)

%tstoolbox/@signal/takens_estimator
%   Syntax:
%     * D2 = takens_estimator2(s, n, range, past)
%
%   Input Arguments:
%     * n - number of randomly chosen reference points (n == -1 means :
%       use all points)
%     * range - maximal relative search radius (relative to attractor
%       size) 0..1
%     * past - number of samples to exclude before and after each
%       reference index
%
%   Takens estimator for correlation dimension
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(4,4)

points = data(s);
[N,dim] = size(points);

try 
	atria = optparams(s, 1);
	names = fieldnames(atria);
catch
	atria = nn_prepare(points, 'euclidian');
	s = setoptparams(s, 1, atria);
end

ref = randref(1,N,n);
D2 = takens_estimator(atria, points, ref, range, past);	

end