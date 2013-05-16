% Calculate negative log-likelihood of the double exponential
% distribution with the rate parameters evaluated at the given data.
%
% Usage:
% nll = exp2like(k1, k2, x)
%
% Inputs
%    k1 - First exponential rate parameter
%    k2 - Second exponential rate parameter
%    x  - Input data
%
% Outputs
%   nll - Negative log-likelihood
%
% (c) Max Little, 2010. If you use this code for your research, please
% cite:
% "Steps and bumps: precision extraction of discrete states of molecular
% machines using physically-based, high-throughput time series analysis"
% Max A. Little et al., 2010, arXiv:1004.1234v1 [q-bio.QM]

function nll = exp2like(k1, k2, x)

error(nargchk(3,3,nargin));

x = x(:);
n = length(x);
nll = -n*log(k1*k2/(k2-k1))-sum(log(exp(-k1*x)-exp(-k2*x)));
