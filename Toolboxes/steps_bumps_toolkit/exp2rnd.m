% Generate double exponentially distributed random values.
%
% Usage:
% x = exp2rnd(k1, k2, N)
%
% Inputs
%    k1 - First exponential rate parameter
%    k2 - Second exponential rate parameter
%    N  - Number of values to generate
%
% Outputs
%    x  - Random values
%
% (c) Max Little, 2010. If you use this code for your research, please
% cite:
% "Steps and bumps: precision extraction of discrete states of molecular
% machines using physically-based, high-throughput time series analysis"
% Max A. Little et al., 2010, arXiv:1004.1234v1 [q-bio.QM]

function x = exp2rnd(k1, k2, N)

error(nargchk(3,3,nargin));

x = exprnd(1/k1,N,1) + exprnd(1/k2,N,1);
