function p = BF_LjungBoxTest(x, maxLag)
% BF_LjungBoxTest   p-values of the Ljung-Box Q-test for autocorrelation at lags 1 to maxLag.
%
% Q(L) = N (N + 2) sum_{k=1}^{L} rho_k^2 / (N - k), with rho_k the sample autocorrelation of
% the mean-removed series, is compared with a chi-squared distribution with L degrees of
% freedom.
%
%---INPUTS:
% x, the series to test (the already-squared series when testing for ARCH effects)
% maxLag, the largest lag (default: 20)
%
%---OUTPUT:
% p, the p-values for lags 1 to maxLag, as a row vector (all NaN for a constant or
%       non-finite series)
%
% Reference: G. M. Ljung and G. E. P. Box, "On a Measure of Lack of Fit in Time Series
% Models", Biometrika 65, 297 (1978). DOI: 10.1093/biomet/65.2.297

% ------------------------------------------------------------------------------
% Copyright (C) 2013-2026, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite the following two papers:
%
% (1) B.D. Fulcher and N.S. Jones, "hctsa: A Computational Framework for Automated
% Time-Series Phenotyping Using Massive Feature Extraction, Cell Systems 5: 527 (2017).
% DOI: 10.1016/j.cels.2017.10.001
%
% (2) B.D. Fulcher, M.A. Little, N.S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013).
% DOI: 10.1098/rsif.2013.0048
%
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program. If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

if nargin < 2 || isempty(maxLag)
    maxLag = 20;
end
x = x(:);
N = length(x);
if maxLag > N - 1
    error('Number of lags must be less than the sample size');
end
p = NaN(1, maxLag);
if ~all(isfinite(x)) || max(x) == min(x)
    return
end

xc = x - mean(x);
den = sum(xc.^2);
rho = zeros(1, maxLag);
for k = 1:maxLag
    rho(k) = sum(xc(k+1:N) .* xc(1:N-k)) / den;
end
Q = N * (N + 2) * cumsum(rho.^2 ./ (N - (1:maxLag)));
p = gammainc(Q / 2, (1:maxLag) / 2, 'upper');

end
