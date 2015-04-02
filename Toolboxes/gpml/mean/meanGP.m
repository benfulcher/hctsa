function A = meanGP(hyp,inf,mean,cov,lik,x,y, hypz,z,i)

% Mean function being the predictive mean of a GP model:
%
% m(z) = posterior mean of another GP at location z as given by
% m(z) = gp(hyp,inf,mean,cov,lik,x,y, z)
%
% The hyperparameters are:
%
% hypz = [ ]
%
% Copyright (c) by Hannes Nickisch, 2014-11-01.
%
% See also MEANFUNCTIONS.M and MEANGPEXACT.M.

if nargin<7, error('GP must be specified.'), end           % check for dimension
if nargin<9, A = '0'; return, end             % report number of hyperparameters

if nargin==9
  A = gp(hyp,inf,mean,cov,lik,x,y, z);                 % evaluate posterior mean
else
  A = zeros(size(x,1),1);                                           % derivative
end
