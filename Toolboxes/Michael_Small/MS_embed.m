% [x,y] or x = MS_embed(z,lags) or MS_embed(z,dim,lag)
% embed z using given lags or dim and lag
% embed(z,dim,lag) == MS_embed(z,[0:lag:lag*(dim-1)])
% negative entries of lags are into future
%
% If return is [x,y], then x is the positive lags and y the negative lags
% Order of rows in x and y the same as sort(lags)
%
% defaults:
%  dim = 3
%  lag = 1
%  lags = [0 1 2]; or [-1 lags] when two outputs and no negative lags
%
% Michael Small
% michael.small@uwa.edu.au, http://school.maths.uwa.edu.au/~small/
% 3/3/2005
% For further details, please see M. Small. Applied Nonlinear Time Series
% Analysis: Applications in Physics, Physiology and Finance. Nonlinear Science
% Series A, vol. 52. World Scientific, 2005. (ISBN 981-256-117-X) and the
% references therein.
% (minor cosmetic changes by Ben Fulcher, 2010)

function [x, y] = MS_embed(z,v,w)

if nargin == 3
  v = (0:w:w*(v-1));
end
if nargin == 1
  v = [0, 1, 2];
end
if nargout == 2 && min(v) >= 0
  v = [-1, v];
end
lags = sort(v);

dim = length(lags);

[c, n] = size(z);
if c ~= 1
  z = z';
  [c, n] = size(z);
end
if c ~= 1
  error('Embed needs a vector as first arg.');
end

if n < lags(dim)
  fprintf(1,'Vector is too small to be embedded with the given lags'); % ++BF changed from error
  x = NaN;
end


w = lags(dim) - lags(1); 		% window
m = n - w; 				% Rows of x
t = (1:m)  + lags(dim); 		% embed times

x = zeros(dim,m);

for i=1:dim
  x(i,:) = z( t  -  lags(i) );
end

if nargout==2
  y = x(v<0,:);
  x = x(v>=0,:);
end

end