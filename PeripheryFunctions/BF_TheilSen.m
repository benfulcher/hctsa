function p = BF_TheilSen(x, y)
% BF_TheilSen   Theil-Sen robust straight-line fit: a closed-form, deterministic alternative to iterative robust regression.
%
% The slope is the median of the slopes of the lines through all pairs of points
% with different x values, and the intercept is the median of y - slope*x. The
% estimate tolerates up to about 29% of points being outliers, and involves no
% iteration, tuning constant or random numbers.
%
% ---INPUTS:
% x, the predictor values (vector)
% y, the response values (vector, same length)
%
% ---OUTPUT:
% p, [slope, intercept], in the same order as polyfit(x, y, 1) (NaN if fewer than
%       two distinct x values).

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

x = x(:);
y = y(:);
n = length(x);

% Slopes of the lines through every pair of points with different x:
[ii, jj] = find(triu(true(n), 1));
dx = x(jj) - x(ii);
good = (dx ~= 0);
if ~any(good)
	p = [NaN, NaN];
	return
end
slope = median((y(jj(good)) - y(ii(good))) ./ dx(good));
p = [slope, median(y - slope*x)];

end
