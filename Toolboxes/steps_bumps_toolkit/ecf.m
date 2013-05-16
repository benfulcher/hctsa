% Calculates empirical characteristic function.
%
% Usage:
% X = ecf(x, f)
%
% Inputs
%    x - Input signal (range 0 to 2pi)
%    f - Fourier analysis frequencies
%
% Outputs
%    X - Empirical characteristic function coefficients
%
% (c) Max Little, 2010. If you use this code for your research, please
% cite:
% "Steps and bumps: precision extraction of discrete states of molecular
% machines using physically-based, high-throughput time series analysis"
% Max A. Little et al., 2010, arXiv:1004.1234v1 [q-bio.QM]

function X = ecf(x, f)

error(nargchk(3,3,nargin));

x = x(:);
f = f(:);

N = length(x);
F = length(f);
X = zeros(F, 1);
for i = 1:F
    X(i) = 1/N*sum(exp(1i*f(i)*x));
end
