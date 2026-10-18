function p = BF_ARCHTest(x, maxLag)
% BF_ARCHTest   p-values of Engle's test for ARCH effects at lags 1 to maxLag.
%
% For each lag L, the squared series (not de-meaned) is regressed on a constant and L of
% its own lags, and the Lagrange-multiplier statistic (N - L) R^2 is compared with a
% chi-squared distribution with L degrees of freedom.
%
%---INPUTS:
% x, the (residual) series, not squared
% maxLag, the largest lag (default: 20)
%
%---OUTPUT:
% p, the p-values for lags 1 to maxLag, as a row vector (all NaN for a constant or
%       non-finite series)
%
% Reference: R. F. Engle, "Autoregressive Conditional Heteroscedasticity with Estimates of the
% Variance of United Kingdom Inflation", Econometrica 50, 987 (1982). DOI: 10.2307/1912773

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
if N - maxLag < 2
    error('The series is too short for the requested number of ARCH lags');
end
p = NaN(1, maxLag);
x2 = x.^2;
if ~all(isfinite(x2)) || max(x2) == min(x2)
    return
end

for L = 1:maxLag
    % Response and lagged regressors over the N - L complete rows
    resp = x2(L+1:N);
    X = ones(N - L, L + 1);
    for k = 1:L
        X(:, k + 1) = x2(L+1-k:N-k);
    end
    beta = X \ resp;
    res = resp - X * beta;
    sst = sum((resp - mean(resp)).^2);
    if sst == 0
        continue
    end
    R2 = 1 - sum(res.^2) / sst;
    p(L) = gammainc((N - L) * R2 / 2, L / 2, 'upper');
end

end
