function [m,dm] = meanGP(hyp,inf,mean,cov,lik,x,y, hypz,z)

% Mean function being the predictive mean of a GP model:
%
% m(z) = posterior mean of another GP at location z as given by
% m(z) = gp(hyp,inf,mean,cov,lik,x,y, z)
%
% The hyperparameters are:
%
% hypz = [ ]
%
% Copyright (c) by Hannes Nickisch, 2016-04-16.
%
% See also meanFunctions.m and meanGPexact.m.

if nargin<7, error('GP must be specified.'), end           % check for dimension
if nargin<9, m = '0'; return, end             % report number of hyperparameters

m = gp(hyp,inf,mean,cov,lik,x,y, z);                   % evaluate posterior mean
dm = @(q) zeros(0,1);                                   % directional derivative