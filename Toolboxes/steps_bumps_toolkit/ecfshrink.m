% Calculates empirical characteristic function, and applies coefficient
% shrinkage, then reconstructs the shrunken density.
%
% Usage:
% [X, Xs, p] = ecfshrink(x, fmax, xk, N, flp)
%
% Inputs
%    x    - Input signal (range 0 to 2pi)
%    fmax - Maximum integer Fourier analysis frequency (non-negative only)
%    xr   - Values of x at which to reconstruct p(x)
%    N    - Number of significant Fourier coefficients to retain;
%           set N=fmax to bypass shrinkage
%    flp  - (Optional) Low pass filtering frequency, set flp=fmax to bypass
%
% Outputs
%    X    - Empirical characteristic function coefficients
%    Xs   - Shrunken characteristic function coefficients
%    p    - Shrunken distribution estimate p(tk) at locations tk
%
% (c) Max Little, 2010. If you use this code for your research, please
% cite:
% "Steps and bumps: precision extraction of discrete states of molecular
% machines using physically-based, high-throughput time series analysis"
% Max A. Little et al., 2010, arXiv:1004.1234v1 [q-bio.QM]

function [X, Xs, p] = ecfshrink(x, fmax, xr, N, flp)

error(nargchk(4,5,nargin));

if (nargin == 4)
    flp = fmax;
end

x = x(:);
xr = xr(:);

X = ecf(x, 0:fmax);
[Xsort,i] = sort(abs(X).^2,'descend');
i = i(1:N);
Xr = zeros(2*fmax+1,1);
Xr(fmax+1) = 1;
Xr(i+fmax) = X(i);
Xr(fmax+2-i) = conj(X(i));

fs = -fmax:fmax;
i = find(abs(fs) > flp);
Xr(i) = 0;
p = ecf2pdf(Xr, fs, xr);
Xs = Xr(fmax+1:end);
