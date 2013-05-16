% Fit double exponential distribution to a set of values, using the method
% of moments.
%
% Usage:
% [k1, k2] = exp2fit(x)
%
% Inputs
%    x  - Input data
%
% Outputs
%    k1 - First exponential rate parameter
%    k2 - Second exponential rate parameter
%
% (c) Max Little, 2010. If you use this code for your research, please
% cite:
% "Steps and bumps: precision extraction of discrete states of molecular
% machines using physically-based, high-throughput time series analysis"
% Max A. Little et al., 2010, arXiv:1004.1234v1 [q-bio.QM]

function [k1, k2] = exp2fit(x)

error(nargchk(1,1,nargin));

x = x(:);

% Method of moments
mu = mean(x);
sg = std(x);
k2 = (mu+sqrt(2*sg^2-mu^2))/(mu^2-sg^2);
k1 = k2/(mu*k2-1);
