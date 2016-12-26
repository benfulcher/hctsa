function [indices, distances] = brute(points, refind, nnr, past)

% [indices, distances] = brute(points, refind, nnr, past)
%
% Brute force implementation of nearest neighbor search
%
% Input arguments :
% 
% points - N by D matrix of N points of dimension D
% refind - integer reference indices
% nnr - number of neighbor to find
% past - number of points to exclude before and after the actual reference index
%
% Output arguments :
%
% indices - R by nnr integer matrix containing the indices of the nearest neighbors, R = length(refind)
% dimension - R by nnr double matrix containing the distances to the nearest neighbors 
%
% Using euclidian norm

narginchk(2,4)

if nnr < 1
	error('At least one neighbor must be requested')
end

[N, dim] =  size(points);
R = length(refind);

indices = zeros(R, nnr);
distances = zeros(R, nnr);

for i=1:R
	actual = refind(i);
	d = euclnorm(points - repmat(points(actual,:), N, 1));
	d(max(1, actual-past):min(N, actual+past)) = inf;
	[d2,ind] = sort(d);
	indices(i,:) = ind(1:nnr)';
	distances(i,:) = d2(1:nnr)';
end


function y = euclnorm(x)

y = sqrt(sum(x.*x, 2));
