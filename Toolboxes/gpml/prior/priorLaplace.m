function [lp,dlp] = priorLaplace(mu,s2,x)

% Univariate Laplacian hyperparameter prior distribution.
% Compute log-likelihood and its derivative or draw a random sample.
% The prior distribution is parameterized as:
%
%   p(x) = exp(-abs(x-mu)/b)/(2*b), where b = sqrt(s2/2),
%
% mu(1x1) is the mean parameter, s2(1x1) is the variance parameter and
% x(1xN) contains query hyperparameters for prior evaluation.
%
% For more help on design of priors, try "help priorDistributions".
%
% Copyright (c) by Roman Garnett and Hannes Nickisch, 2014-09-09.
%
% See also priorDistributions.m.

if nargin<3                                                    % return a sample
  u = rand-1/2; lp = mu - sqrt(s2/2)*sign(u)*log(1-2*abs(u)); return
end

b = sqrt(s2/2);
lp  = -abs(x-mu)/b - log(2*b);
dlp = -sign(x-mu)/b;