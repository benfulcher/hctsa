function A = meanOne(hyp, x, i)

% One mean function. The mean function does not have any parameters.
%
% m(x) = 1
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2010-08-04.
%
% See also MEANFUNCTIONS.M.

if nargin<2, A = '0'; return; end             % report number of hyperparameters
if nargin==2
  A = ones(size(x,1),1);                                         % evaluate mean
else
  A = zeros(size(x,1),1);                                           % derivative
end