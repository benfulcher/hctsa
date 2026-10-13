function f = BF_ExpFit(x, y, withOffset, maxRate)
% BF_ExpFit   Global least-squares fit of an exponential, a*exp(b*x) + c, by variable projection.
%
% For a given rate b the best amplitude a (and offset c) follow from linear least
% squares in closed form, so the fit reduces to a search over the single rate b.
% This is done on a fixed grid of rates (always including b = 0, a constant), then
% refined by 40 steps of golden-section search around the best grid point. The
% result is the global optimum within the allowed range of rates (rather than the
% local optimum an iterative nonlinear fit stops at), and it needs no starting
% point. When the curve is close to straight, the rate settles near zero while a
% and c become large and opposite in sign, so a and c are poorly determined
% (whereas b, r2 and rmse are not). If the best rate lies at the edge of the allowed
% range (the data call for a step rather than an exponential), b is returned as that
% limiting value, i.e. "at least this steep".
%
% ---INPUTS:
% x, the predictor values (vector)
% y, the data to fit (vector, same length)
% withOffset, true to fit a*exp(b*x) + c, false to fit a*exp(b*x) (default: true)
% maxRate, the largest allowed |b| in units of 1/range(x): rates are searched
%       between -maxRate/range(x) and maxRate/range(x) (default: 20, so the fitted
%       curve changes by at most a factor of exp(20) across the data)
%
% ---OUTPUT: a structure with fields
% a, b, c: the parameters of the fit (c = 0 if there is no offset)
% r2: the coefficient of determination, 1 - SSE/SST, which lies in [0,1] because
%       the search includes the constant model b = 0
% adjr2: adjusted R^2, 1 - (1 - r2)*(n - 1)/(n - p), for p = 3 (or 2) parameters
% rmse: root-mean-square error of the fit, sqrt(SSE/(n - p))
% All are NaN if y is constant or not finite, or if there are not more data points
% than parameters.

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

if nargin < 3 || isempty(withOffset)
	withOffset = true;
end
if nargin < 4 || isempty(maxRate)
	maxRate = 20;
end

x = x(:);
y = y(:);
n = length(y);
numParams = 2 + withOffset;
yc = y - mean(y);
SST = sum(yc.^2);
if n <= numParams || ~(SST > 0) || ~all(isfinite(y)) || range(x) == 0
	f = struct('a', NaN, 'b', NaN, 'c', NaN, 'r2', NaN, 'adjr2', NaN, 'rmse', NaN);
	return
end

% Search the rate on a grid, then refine around the best grid point:
bMax = maxRate / range(x);
bGrid = linspace(-bMax, bMax, 401);
[~, k] = min(SSEofRate(bGrid, x, y, withOffset));
lo = bGrid(max(k - 1, 1));
hi = bGrid(min(k + 1, length(bGrid)));
phi = (sqrt(5) - 1) / 2;
for i = 1:40 % golden-section search
	b1 = hi - phi*(hi - lo);
	b2 = lo + phi*(hi - lo);
	if SSEofRate(b1, x, y, withOffset) <= SSEofRate(b2, x, y, withOffset)
		hi = b2;
	else
		lo = b1;
	end
end
b = (lo + hi) / 2;

% Closed-form amplitude and offset at the chosen rate:
e = exp(b*x);
if withOffset
	a = sum((e - mean(e)).*yc) / sum((e - mean(e)).^2);
	c = mean(y) - a*mean(e);
else
	a = sum(e.*y) / sum(e.^2);
	c = 0;
end
SSE = sum((y - a*e - c).^2);

f.a = a;
f.b = b;
f.c = c;
f.r2 = min(max(1 - SSE/SST, 0), 1); % (the constant model is within the search: only rounding can leave [0,1])
f.adjr2 = 1 - (1 - f.r2)*(n - 1)/(n - numParams);
f.rmse = sqrt(SSE/(n - numParams));

end

% ------------------------------------------------------------------------------
function sse = SSEofRate(bs, x, y, withOffset)
	% Sum of squared errors of the best-fitting a (and c) at each rate in bs
	e = exp(x*bs); % each column is exp(b*x) for one rate
	if withOffset
		ec = e - mean(e, 1);
		den = sum(ec.^2, 1);
		sse = sum((y - mean(y)).^2) - (ec'*(y - mean(y)))'.^2 ./ den;
	else
		den = sum(e.^2, 1);
		sse = sum(y.^2) - (e'*y)'.^2 ./ den;
	end
	sse(~(den > 0) | ~isfinite(sse)) = Inf; % (b = 0 has no offset-model solution)
end
