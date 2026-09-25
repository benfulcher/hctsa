function out = CO_Embed2_Shapes(y, tau, shape, r, theilerWin)
% CO_Embed2_Shapes Shape-based statistics in a 2-d embedding space
%
% Takes a shape and places it on each point in the two-dimensional time-delay
% embedding space sequentially. This function counts the points inside this shape
% as a function of time, and returns statistics on this extracted time series.
%
% ---INPUTS:
% y, the input time-series as a (z-scored) column vector
% tau, the time-delay
% shape, has to be 'circle' for now...
% r, the radius of the circle
% theilerWin, Theiler window: points closer in time than this are not counted
%       as neighbors, since they are close in the embedding only because
%       successive values are correlated ({'ac', k} for k times the first
%       zero-crossing of the autocorrelation function, or a number of samples;
%       see BF_TheilerWindow; default: {'ac', 1}; 0 counts every other point)
%
% ---OUTPUTS:
% The constructed time series of the number of nearby points, and
% include the autocorrelation, spread (std, IQR), median, mode, histogram
% entropy, and stationarity over fifths of the time series.
%
% ---NOTES:
% `max` was dropped 2026-08-11: redundancy-checked against `mean`/`std` on
% Bonn EEG (500 series) and Empirical1000 (1000 series), |r|>=0.9 with both
% on both datasets, in both registered radii (r=0.1 and r=1).

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

doPlot = 0; % plot results for debugging

% ------------------------------------------------------------------------------
%% Check inputs, set defaults:
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(tau)
	fprintf(1, ['Setting tau as the first zero crossing ' ...
				'of the autocorrelation function.\n']);
	tau = 'tau';
end
if nargin < 3 || isempty(shape)
	shape = 'circle'; % use a circle by default
	fprintf(1, 'Using a circle.\n');
end
if nargin < 4 || isempty(r)
	r = 1; % default radius of 1
end
if nargin < 5 || isempty(theilerWin)
	theilerWin = {'ac', 1};
end
theilerWin = BF_TheilerWindow(y, theilerWin);
if isnan(theilerWin) % the autocorrelation function never crosses zero
	out = NaN; return
end

% Can set time lag equal to first zero crossing of the autocorrelation function with the 'tau' input
if strcmp(tau, 'tau'),
	tau = CO_FirstCrossing(y, 'ac', 0, 'discrete');
	if isnan(tau)
		out = NaN; return
	end
	% Cannot set the time delay greater than 10% the length of the time series
	if tau > length(y) / 10
		tau = floor(length(y) / 10);
	end
end

% Ensure y is a column vector:
if size(y, 2) > size(y, 1);
	y = y';
end

% ------------------------------------------------------------------------------
%% Create the recurrence space, populated by points m
% ------------------------------------------------------------------------------
m = [y(1:end - tau), y(1 + tau:end)];
N = length(m);

if doPlot % plot the recurrence space:
	plot(m(:, 1), m(:, 2), '.');
end

%% Start the analysis

counts = zeros(N, 1); % stores the counts for points within circle

switch shape
	case 'circle'
		% Puts a circle around each point in the embedding space in turn
		% counts how many points are inside this shape, looks at the time series thus formed

		for i = 1:N % across all points in the time series

			% m - m(i,:) relies on implicit broadcasting to subtract the
			% row m(i,:) from every row of m -- identical result to
			% m - ones(N,1)*m(i,:), but without allocating/multiplying by an
			% explicit ones(N,1) matrix on every iteration. (A single fully
			% vectorized N x N pairwise-distance computation across all i
			% would avoid the loop entirely, but would need an N x N x 2
			% intermediate array -- risking excessive memory use for long
			% time series -- so the per-i loop is kept.)
			m_c = m - m(i, :); % points wrt current point i
			m_c_d = sum(m_c.^2, 2); % Euclidean distances from point i

			% Fraction of the points outside i's Theiler window (which includes
			% i itself) that are enclosed in a circle of radius r:
			isOther = abs((1:N)' - i) > theilerWin;
			counts(i) = sum(m_c_d(isOther) <= r^2) / sum(isOther);
		end

	otherwise
		error('Unknown shape ''%s''', shape)
end
% ------------------------------------------------------------------------------
% No point has any (non-Theiler-excluded) neighbor within r
% ------------------------------------------------------------------------------
% A zero density everywhere is a valid (if extreme) outcome: the location/spread
% statistics are zero, and statistics of the shape of the (constant) counts are
% undefined (NaN)
if all(counts == 0)
	out = struct('ac1', NaN, 'ac2', NaN, 'ac3', NaN, 'tau', NaN, 'std', 0, ...
				'median', 0, 'mean', 0, 'iqr', 0, 'iqronrange', NaN, ...
				'mode_val', 1, 'mode', 0, 'hist_ent', 0, 'statav5_m', NaN, 'statav5_s', NaN);
	return
end

% (counts are the fraction of the embedding within r of each point, i.e. the
% pointwise correlation sum, which is intensive and bounded in [0,1]: raw counts
% would grow in direct proportion to N)

% ------------------------------------------------------------------------------
% Return basic statistics on the counts
% ------------------------------------------------------------------------------
out.ac1 = CO_AutoCorr(counts, 1, 'Fourier');
out.ac2 = CO_AutoCorr(counts, 2, 'Fourier');
out.ac3 = CO_AutoCorr(counts, 3, 'Fourier');
out.tau = CO_FirstCrossing(counts, 'ac', 0, 'continuous');
out.std = std(counts);
out.median = median(counts);
out.mean = mean(counts);
out.iqr = iqr(counts);
out.iqronrange = out.iqr / range(counts);

% --- Distribution
% Using the sqrt binning method:
[binP, binEdges] = histcounts(counts, 'BinMethod', 'sqrt', 'Normalization', 'probability');
binCentres = mean([binEdges(1:end - 1); binEdges(2:end)]);
[out.mode_val, mix] = max(binP);
out.mode = binCentres(mix);
% --- histogram entropy:
out.hist_ent = sum(binP(binP > 0) .* log(binP(binP > 0)));

if doPlot
	bar(binCentres, binP);
end

% --- Stationarity measure for fifths of the time series
afifth = floor(N / 5);
buffer_m = zeros(afifth, 5); % stores a fifth of the time series (embedding vector) in each entry
for i = 1:5
	buffer_m(:, i) = counts(afifth * (i - 1) + 1:afifth * i);
end
out.statav5_m = std(mean(buffer_m)) / std(counts);
out.statav5_s = std(std(buffer_m)) / std(counts);

end
