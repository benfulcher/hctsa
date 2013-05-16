% Reconstruct probability density from characteristic function coefficients.
%
% Usage:
% p = ecf2pdf(X, f, x)
%
% Inputs
%    X - Characteristic function coefficients
%    f - Frequencies associated with each X
%    x - Values of x at which to reconstruct p(x)
%
% Outputs
%    p - Estimated values of p(x) at each given x
%
% (c) Max Little, 2010. If you use this code for your research, please
% cite:
% "Steps and bumps: precision extraction of discrete states of molecular
% machines using physically-based, high-throughput time series analysis"
% Max A. Little et al., 2010, arXiv:1004.1234v1 [q-bio.QM]

function p = ecf2pdf(X, f, x)

error(nargchk(3,3,nargin));

x = x(:);

N = length(x);
p = zeros(N, 1);
for i = 1:N
    p(i) = real(1/N*sum(exp(-1i*x(i)*f).*conj(X)'));
end
