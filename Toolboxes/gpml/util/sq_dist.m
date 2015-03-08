% sq_dist - a function to compute a matrix of all pairwise squared distances
% between two sets of vectors, stored in the columns of the two matrices, a
% (of size D by n) and b (of size D by m). If only a single argument is given
% or the second matrix is empty, the missing matrix is taken to be identical
% to the first.
%
% Usage: C = sq_dist(a, b)
%    or: C = sq_dist(a)  or equiv.: C = sq_dist(a, [])
%
% Where a is of size Dxn, b is of size Dxm (or empty), C is of size nxm.
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2010-12-13.

function C = sq_dist(a, b)

if nargin<1  || nargin>3 || nargout>1, error('Wrong number of arguments.'); end
bsx = exist('bsxfun','builtin');      % since Matlab R2007a 7.4.0 and Octave 3.0
if ~bsx, bsx = exist('bsxfun'); end      % bsxfun is not yet "builtin" in Octave
[D, n] = size(a);

% Computation of a^2 - 2*a*b + b^2 is less stable than (a-b)^2 because numerical
% precision can be lost when both a and b have very large absolute value and the
% same sign. For that reason, we subtract the mean from the data beforehand to
% stabilise the computations. This is OK because the squared error is
% independent of the mean.
if nargin==1                                                     % subtract mean
  mu = mean(a,2);
  if bsx
    a = bsxfun(@minus,a,mu);
  else
    a = a - repmat(mu,1,size(a,2));  
  end
  b = a; m = n;
else
  [d, m] = size(b);
  if d ~= D, error('Error: column lengths must agree.'); end
  mu = (m/(n+m))*mean(b,2) + (n/(n+m))*mean(a,2);
  if bsx
    a = bsxfun(@minus,a,mu); b = bsxfun(@minus,b,mu);
  else
    a = a - repmat(mu,1,n);  b = b - repmat(mu,1,m);
  end
end

if bsx                                               % compute squared distances
  C = bsxfun(@plus,sum(a.*a,1)',bsxfun(@minus,sum(b.*b,1),2*a'*b));
else
  C = repmat(sum(a.*a,1)',1,m) + repmat(sum(b.*b,1),n,1) - 2*a'*b;
end
C = max(C,0);          % numerical noise can cause C to negative i.e. C > -1e-14
