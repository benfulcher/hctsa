% Calculate mean and standard deviation of a double exponential
% distribution.
%
% Usage:
% [m, s] = exp2stat(k1, k2)
%
% Inputs
%    k1 - First exponential rate parameter
%    k2 - Second exponential rate parameter
%
% Outputs
%    m  - Mean
%    s  - Standard deviation
%
% (c) Max Little, 2010. If you use this code for your research, please
% cite:
% "Steps and bumps: precision extraction of discrete states of molecular
% machines using physically-based, high-throughput time series analysis"
% Max A. Little et al., 2010, arXiv:1004.1234v1 [q-bio.QM]

function [m, s] = exp2stat(k1, k2)

error(nargchk(2,2,nargin));

m = (k1+k2)/(k1*k2);
s = sqrt((k1^2+k2^2)/(k1^2*k2^2));
