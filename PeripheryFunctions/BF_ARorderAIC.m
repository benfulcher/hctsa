function p = BF_ARorderAIC(y, maxLag)
% BF_ARorderAIC   The order of a zero-mean AR model selected by Akaike's information criterion.
%
% AR models of orders 0, 1, ..., L (L = min(20, floor(N/10)), at least 1, unless given)
% are fitted by ordinary least squares on the same samples, y(L+1:N), and the order
% minimizing AIC(p) = log(s2(p)) + 2p/n is returned, where n = N - L and s2(p) is the mean
% squared residual (the mean square of the data for p = 0). The variance s2(p) is floored
% at 1e-10 times the mean square of the data, so that orders whose residuals are at the
% level of rounding error (exactly predictable series) tie, and the penalty then picks
% the lowest of them.
%
%---NOTES:
% H. Akaike, "A new look at the statistical model identification", IEEE Trans. Autom.
% Control 19, 716 (1974).

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

y = y(:);
N = length(y);
if nargin < 2 || isempty(maxLag)
    maxLag = max(min(20, floor(N/10)), 1);
end
L = maxLag;
if N - L < L + 2
    p = 0; return
end
target = y(L+1:N);
n = length(target);
H = hankel(y(1:L+1), y(L+1:N))'; % row j = y(j:j+L)
lags = H(:, L:-1:1); % column k = y(t-k)
floorS2 = 1e-10 * mean(target.^2);
aic = zeros(L + 1, 1);
for k = 0:L
    if k == 0
        res = target;
    else
        X = lags(:, 1:k);
        res = target - X * (BF_Pinv(X) * target);
    end
    s2 = max(mean(res.^2), floorS2);
    aic(k + 1) = log(s2) + 2 * k / n;
end
[~, ix] = min(aic);
p = ix - 1;

end
