function z = BF_RunsZ(y)
% BF_RunsZ   Signed z-statistic of a runs test for randomness about the median.
%
% Splits the values of y into those above the median and those at or below it,
% counts the runs (maximal stretches of consecutive values on the same side) and
% standardizes the count with its exact mean and variance under the null
% hypothesis that the values occur in random order (Wald and Wolfowitz):
%
%   z = (R - mu) / sigma,  mu = 1 + 2*n1*n2/n,
%   sigma^2 = 2*n1*n2*(2*n1*n2 - n) / (n^2 * (n - 1)),
%
% where R is the number of runs, n1 and n2 are the numbers of values above and
% at or below the median, and n = n1 + n2.
%
% The statistic is returned instead of a p-value because the p-value of a series
% with strong serial structure underflows to (numerically) zero, and then only
% its order of magnitude carries information. The z-statistic keeps ordering such
% series, and its sign gives the direction of the departure from randomness:
% z < 0 (too few runs) indicates positive serial dependence (slowly varying
% series, trends), z > 0 (too many runs) indicates negative serial dependence
% (alternation). Values at the median are assigned to the lower group (not
% discarded), so a series with many tied values is still tested. If no value
% lies above the median (more than half of the values equal the maximum), the
% groups are instead the values at the median and those below it.
%
% ---INPUTS:
% y, a vector (NaN values are ignored)
%
% ---OUTPUT:
% z, the standardized number of runs; approximately standard normal for a
%    random ordering of the values. NaN if no runs test is possible: for a
%    constant series, or when the null distribution of the number of runs has
%    zero variance (two values, one on each side of the median).

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

y = y(~isnan(y(:)));
m = median(y);
isUp = y > m; % above the median (values at the median count as below)
if ~any(isUp)
	isUp = y >= m; % more than half the values are at the maximum: split at the median instead
end
n1 = sum(isUp); % number of values in the upper group
n2 = length(y) - n1; % number of values in the lower group
if n1 == 0 || n2 == 0
	z = NaN; % no values on one side: a runs test is not possible
	return
end
n = n1 + n2;

R = 1 + sum(isUp(2:end) ~= isUp(1:end - 1)); % number of runs
mu = 1 + 2*n1*n2/n; % expected number of runs
v = 2*n1*n2*(2*n1*n2 - n) / (n^2 * (n - 1)); % variance of the number of runs
if v == 0
	z = NaN; % only n1 = n2 = 1: the number of runs is fixed
else
	z = (R - mu) / sqrt(v);
end

end
