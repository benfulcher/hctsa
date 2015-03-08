function A = meanNN(c,m, hyp, x, i)

% Nearest neighbor mean function. The mean function is parameterized as:
%
% m(z) = m_j, j = arg min_i d(ci,x) where d is the Euclidean distance and ci is
%                                   the ith cluster center.
%
% The hyperparameters are:
%
% hyp = [ ]
%
% Copyright (c) by Hannes Nickisch, 2014-12-04.
%
% See also MEANFUNCTIONS.M.

if nargin<4, A = '0'; return; end             % report number of hyperparameters 
if numel(hyp)~=0, error('No hyperparameters needed for this model.'), end
if nargin==4                                                     % evaluate mean
  [junk,j] = min(sq_dist(c',x'));
  A = m(j); A = A(:);
else
  A = zeros(size(x,1),1);                                           % derivative
end