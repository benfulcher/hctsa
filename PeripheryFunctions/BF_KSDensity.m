function [f, xi, h] = BF_KSDensity(y, xi, h)
% BF_KSDensity   Gaussian kernel density estimate with an explicit bandwidth.
%
% Evaluates the kernel-smoothed density of the data, f(x) = mean_i N(x; y_i, h^2),
% exactly (no binning and no truncation of the kernel), with the bandwidth written
% out as a formula instead of left to the defaults of ksdensity. The default
% bandwidth is Silverman's normal-reference rule with a robust scale estimate:
%       h = s * (4/(3*n))^(1/5),   s = median(|y - median(y)|)/0.6745,
% where s estimates the standard deviation from the median absolute deviation (this
% is also the default of ksdensity). If the median absolute deviation is zero (more
% than half the values are equal), s is the standard deviation of y instead; if that
% too is zero (a constant series), there is no scale to smooth over and h, xi and f
% are NaN.
%
% ---INPUTS:
% y, the data vector (NaNs are ignored)
% xi, the points at which to evaluate the density (default: 100 equally spaced
%       points from min(y) - 3*h to max(y) + 3*h)
% h, the kernel bandwidth (standard deviation of the Gaussian kernel; default: the
%       rule above)
%
% ---OUTPUTS:
% f, the density estimate at xi (integrates to 1 over the real line)
% xi, the points at which it is evaluated (a row vector)
% h, the bandwidth used

% ------------------------------------------------------------------------------
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

y = y(~isnan(y));
y = y(:);
n = length(y);

% ------------------------------------------------------------------------------
% Bandwidth
% ------------------------------------------------------------------------------
if nargin < 3 || isempty(h)
	s = median(abs(y - median(y))) / 0.6745; % robust estimate of the standard deviation
	if s <= 0
		s = std(y);
	end
	if s <= 0
		h = NaN; % constant data: no scale
	else
		h = s * (4 / (3 * n))^(1/5);
	end
end

% ------------------------------------------------------------------------------
% Evaluation points
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(xi)
	xi = linspace(min(y) - 3 * h, max(y) + 3 * h, 100);
end
xi = reshape(xi, 1, []);

% ------------------------------------------------------------------------------
% Sum the Gaussian kernels (in blocks of evaluation points, to limit memory)
% ------------------------------------------------------------------------------
f = zeros(size(xi));
blockSize = max(1, floor(2e6 / n));
for i = 1:blockSize:numel(xi)
	ix = i:min(i + blockSize - 1, numel(xi));
	f(ix) = sum(exp(-0.5 * ((xi(ix) - y) / h).^2), 1) / (n * h * sqrt(2 * pi));
end

end
