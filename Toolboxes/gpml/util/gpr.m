function varargout = gpr(hyper, cov, x, y, xs)

% gpr - Gaussian process regression, with a named covariance function. Two
% modes are possible: training and prediction: if no test data are given, the
% function returns minus the log likelihood and its partial derivatives with
% respect to the hyperparameters; this mode is used to fit the hyperparameters.
% If test data are given, then (marginal) Gaussian predictions are computed,
% whose mean and variance are returned. Note that in cases where the covariance
% function has noise contributions, the variance returned in s2 is for noisy
% test targets; if you want the variance of the noise-free latent function, you
% must substract the noise variance.
%
% usage: [nlZ dnlZ] = gpr(hyp, cov, x, y)
%    or: [mu s2]    = gpr(hyp, cov, x, y, xs)
%
% where:
%
%   hyp      is a (column) vector of log hyperparameters
%   cov      is the covariance function
%   x        is a n by D matrix of training inputs
%   y        is a (column) vector (of size n) of targets
%   xs       is a ns by D matrix of test inputs
%   nlZ      is the returned value of the negative log marginal likelihood
%   dnlZ     is a (column) vector of partial derivatives of the negative
%                 log marginal likelihood wrt each log hyperparameter
%   mu       is a (column) vector (of size nn) of prediced means
%   s2       is a (column) vector (of size nn) of predicted variances
%
% For more help on covariance functions, see "help covFunctions".
%
% Copyright (c) 2010 Carl Edward Rasmussen and Hannes Nickisch 2010-06-18.

err = 'we need cov = {''covSum'', {cov1, ..,''covNoise'', .., covP} }; to map to gp.m';
if ~strcmp(cov{1},'covSum'), error(err), end
id = 0; nhyp = zeros(length(cov{2}),1);
for i=1:length(cov{2})
  if strcmp(cov{2}{i},'covNoise'), id = i; break, end
  nhyp(i) = eval(feval(cov{2}{i}));
end
if id==0, error(err), else nhyp = sum(nhyp(1:id-1)); end
cov{2} = cov{2}([1:id-1,id+1:end]);
hyp.lik = hyper(nhyp+1); hyp.cov = hyper([1:nhyp,nhyp+2:end]);

% Note, this function is just a wrapper provided for backward compatibility,
% the functionality is now provided by the more general gp function.
lik = @likGauss; inf = @infExact; mean = @meanZero;

varargout = cell(nargout, 1);    % allocate the right number of output arguments
if nargin==4
  [varargout{:}] = gp(hyp,inf,mean,cov,lik,x,y);
  if nargout>1, varargout{2} = [varargout{2}.lik; varargout{2}.cov]; end
else
  [varargout{:}] = gp(hyp,inf,mean,cov,lik,x,y,xs);
end
