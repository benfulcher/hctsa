function [a, b] = DK_defit(delta,epsilon,maxdelta)
% [a,b] = defit(delta,epsilon,maxdelta)
% linear fitting routine for delta-epsilon
% delta   -- distances between pre-images: output by delta-epsilon
% epsilon -- distances between images: output by delta-epsilon
% maxdelta-- optional - largest delta to consider.
% Copyright (c) 1996 by D. Kaplan, All Rights Reserved

if nargin < 3
  maxdelta = 10e100;
end

% get rid of points whose pre-images are too far apart.
foo = (delta < maxdelta);
delta = delta(foo);
epsilon = epsilon(foo);

% Sort the pairs into ascending order on delta
[x, i] = sort(delta);
y = epsilon(i);

s = length(x);
sx = sum(x);
sy = sum(y);
sxx = sum( x.*x);
sxy = sum( x.*y);
del = s*sxx - sx*sx;
a = (sxx*sy - sx*sxy)/del;
b = (s*sxy - sx*sy)/del;

end