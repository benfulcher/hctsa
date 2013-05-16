% Calculates double exponential density function at given values.
%
% Usage:
% p = exp2pdf(x, k1, k2)
%
% Inputs
%    x  - Input data
%    k1 - First exponential rate parameter
%    k2 - Second exponential rate parameter
%
% Outputs
%    p  - Probability density function evaluated at each x
%
% (c) Max Little, 2010. If you use this code for your research, please
% cite:
% "Steps and bumps: precision extraction of discrete states of molecular
% machines using physically-based, high-throughput time series analysis"
% Max A. Little et al., 2010, arXiv:1004.1234v1 [q-bio.QM]

function p = exp2pdf(x, k1, k2)

error(nargchk(3,3,nargin));

x = x(:);

p = k1*k2/(k2-k1)*(exp(-k1*x)-exp(-k2*x));
