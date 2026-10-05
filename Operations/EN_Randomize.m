function out = EN_Randomize(y, randomizeHow, randomSeed)
% EN_Randomize   How properties of the series change as it is progressively randomized.
%
% Randomizes a copy of the input (z-scored) series one point at a time, according
% to a randomization procedure, repeated 2N times for a series of length N, and
% compares statistics of the randomized copy with the original at 21 checkpoints:
% at the start and after every N/10 steps.
%
% ---INPUTS:
% y, the input (z-scored) time series
% randomizeHow, what one step of randomization does:
%       'statdist': overwrites a random element of the series with a randomly
%                   chosen element of the original series
%       'dyndist': overwrites a random element of the series with another random
%                  element of the current, partially randomized series
%       'permute': swaps two randomly chosen elements of the series, so that the
%                  distribution of values never changes and only the temporal
%                  properties do
%       Default: 'statdist'.
% randomSeed, the seed of the random choices (see BF_RandomSeed; they come from the
%       portable generator BF_Random, so results are reproducible across languages)
%
% ---OUTPUTS: for each of ten statistics measured at each checkpoint, six or seven
% fields describing its trajectory over the 21 checkpoints. The statistics are:
% xcn1, xc1: the cross-correlation of the original and randomized series at lags
%       -1 and +1
% d1: the distance between the original and randomized series, norm(y - y_rand) / N
% ac1, ac2, ac3, ac4: the autocorrelation of the randomized series at lags 1 to 4
% permen3_1: the normalized permutation entropy PermEn(3,1) of the randomized series
% statav5: StatAv with 5 segments (the standard deviation of the segment means)
% swss5_1: the standard deviation across 5 non-overlapping windows of the local
%       standard deviation, relative to the full-series standard deviation
% The fields are named by joining a statistic's name to a suffix:
% xcn1, xc1, ac1, ac2, ac3, ac4 (fits of a * exp(b * k), k the checkpoint number
% 1..21) have the suffixes fexpa, fexpb (the parameters a and b), fexpr2 (R^2),
% fexprmse (root-mean-square error), diff and hp:
%       xcn1fexpa, xcn1fexpb, xcn1fexpr2, xcn1fexprmse, xcn1diff, xcn1hp, xc1fexpa,
%       xc1fexpb, xc1fexpr2, xc1fexprmse, xc1diff, xc1hp, ac1fexpa, ac1fexpb,
%       ac1fexpr2, ac1fexprmse, ac1diff, ac1hp, ac2fexpa, ac2fexpb, ac2fexpr2,
%       ac2fexprmse, ac2diff, ac2hp, ac3fexpa, ac3fexpb, ac3fexpr2, ac3fexprmse,
%       ac3diff, ac3hp, ac4fexpa, ac4fexpb, ac4fexpr2, ac4fexprmse, ac4diff, ac4hp
% d1, permen3_1, statav5, swss5_1 (fits of a * exp(b * k) + c) have the same
% suffixes plus fexpc (the offset c):
%       d1fexpa, d1fexpb, d1fexpc, d1fexpr2, d1fexprmse, d1diff, d1hp,
%       permen3_1fexpa, permen3_1fexpb, permen3_1fexpc, permen3_1fexpr2,
%       permen3_1fexprmse, permen3_1diff, permen3_1hp, statav5fexpa, statav5fexpb,
%       statav5fexpc, statav5fexpr2, statav5fexprmse, statav5diff, statav5hp,
%       swss5_1fexpa, swss5_1fexpb, swss5_1fexpc, swss5_1fexpr2, swss5_1fexprmse,
%       swss5_1diff, swss5_1hp
% In all cases diff is the absolute change |s_end - s_start| of the statistic between
% the first and last checkpoints and hp is the number of the first checkpoint at which
% the statistic passes halfway between its start and end values.
%
% ---NOTES:
% Requires the Curve Fitting Toolbox.
% diff is an absolute change, not a change relative to the starting value, because
% the starting value (e.g., the autocorrelation of the original series at lag 2) can be
% near 0, where a relative change is unstable. All the statistics are on a bounded,
% dimensionless scale (correlations, a normalized entropy, and standard deviations
% of the z-scored series), so the absolute change is comparable across series.

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

% ------------------------------------------------------------------------------
% Check toolboxes, and a z-scored time series:
% ------------------------------------------------------------------------------

% Check a curve-fitting toolbox license is available:
BF_CheckToolbox('curve_fitting_toolbox');

% Check input time series is z-scored:
if ~BF_iszscored(y)
	warning('The input time series should be z-scored for EN_Randomize.')
end

% ------------------------------------------------------------------------------
%% Check inputs:
% ------------------------------------------------------------------------------
% randomizeHow, how to do the randomization:
if nargin < 2 || isempty(randomizeHow)
	randomizeHow = 'statdist'; % use statdist by default
end

% randomSeed: how to treat the randomization
if nargin < 3
	randomSeed = []; % default
end
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
% Preliminaries
% ------------------------------------------------------------------------------

doPlot = false; % Don't plot to screen by default:
N = length(y); % length of the time series

% Set up the points through the randomization process at which the
% calculation of stats will occur:
randp_max = 2; % time series has been randomized to double its length
rand_inc = 0.1; % this proportion of the time series has been randomized between calculations
numCalcs = randp_max / rand_inc; % number of calculations required
calc_ints = floor(randp_max * N / numCalcs);
if calc_ints == 0, calc_ints = 1; end % round up for short time series
calc_pts = (0:calc_ints:randp_max * N);
if calc_pts(end) ~= randp_max * N;
	calc_pts = [calc_pts, randp_max * N];
end
numCalcs = length(calc_pts); % some rounding issues inevitable

statNames = {'xcn1', 'xc1', 'd1', 'ac1', 'ac2', 'ac3', 'ac4', 'permen3_1', 'statav5', 'swss5_1'};
numStats = length(statNames);
stats = zeros(numCalcs, numStats); % record a stat at each randomization increment

y_rand = y; % this vector will be randomized

stats(1, :) = CalculateStats(y, y_rand); % initial condition: apply on itself

% The random choices for every step, reproducible from the seed (portable generator
% BF_Random): two uniform draws per step, as indices uniform on 1..N
randIdx = 1 + floor(N * reshape(BF_Random(2 * randp_max * N, BF_RandomSeed(randomSeed)), 2, []));

% -------------------------------------------------------------------------------
% Do the randomization
% -------------------------------------------------------------------------------
% fprintf(1,'%u/%u calculation points',numCalcs,N*randp_max)

for i = 1:N * randp_max
	switch randomizeHow
		case 'statdist'
			% randomize by substituting a random element of the time series by
			% a random element from the static original time series distribution
			y_rand(randIdx(1, i)) = y(randIdx(2, i));

		case 'dyndist'
			% randomize by substituting a random element of the time series
			% by a random element of the current, already partially randomized,
			% time series
			y_rand(randIdx(1, i)) = y_rand(randIdx(2, i));

		case 'permute'
			% randomize by swapping elements of the time series so that
			% the distribution remains static; only temporal properties will change
			randis = randIdx(:, i);
			tmp = y_rand(randis(1));
			y_rand(randis(1)) = y_rand(randis(2));
			y_rand(randis(2)) = tmp;

		otherwise
			error('Unknown randomization method ''%s''.', randomizeHow);
	end

	if any(calc_pts == i)
		stats(calc_pts == i, :) = CalculateStats(y, y_rand);
	end

end
% fprintf(1,'Randomization took %s',BF_TheTime(toc(randTimer)));

if doPlot
	f = figure('color', 'w'); box('on');
	plot(stats, '.-');
end

% ------------------------------------------------------------------------------
%% Fit exponentials to outputs:
% ------------------------------------------------------------------------------
r = (1:size(stats, 1))'; % gives an 'x-axis' for fitting

% 1) xcn1: cross correlation at lag of -1
% 2) xc1: cross correlation at lag 1
% 3) d1: norm of differences between original and randomized time series
% 4) ac1
% 5) ac2
% 6) ac3
% 7) ac4
% 8) normalized permutation entropy, PermEn(3,1)
% 9) statav5
% 10) swss5_1

out = struct();
for i = 1:length(statNames)
	% Exponential fits:
	switch statNames{i}
		case {'xcn1', 'xc1'}
			startPoint = [stats(1, i), -0.1];
			[c, gof] = f_fix_exp(r, stats(:, i), startPoint, 0);
		case {'ac1', 'ac2', 'ac3'}
			startPoint = [stats(1, i), -0.2];
			[c, gof] = f_fix_exp(r, stats(:, i), startPoint, 0);
		case 'ac4'
			startPoint = [stats(1, i), -0.4];
			[c, gof] = f_fix_exp(r, stats(:, i), startPoint, 0);
		case {'d1', 'permen3_1'}
			startPoint = [-stats(end, i), -0.2, stats(end, i)];
			[c, gof] = f_fix_exp(r, stats(:, i), startPoint, 1);
		case {'statav5', 'swss5_1'}
			startPoint = [-stats(end, i), -0.1, stats(end, i)];
			[c, gof] = f_fix_exp(r, stats(:, i), startPoint, 1);
	end
	out = assignExpStats(out, c, gof, statNames{i});

	% Extra statistics:
	out = assignExtraStats(out, stats(:, i), statNames{i});
end

% ------------------------------------------------------------------------------
function out = CalculateStats(y, y_rand)
	% Calculate statistics comparing a time series, y, and a randomized
	% version of it, y_rand

	% Cross Correlation to original signal
	xc = xcorr(y, y_rand, 1, 'coeff');
	xcn1 = xc(1);
	xc1 = xc(3);

	% Norm of differences between original and randomized signals
	d1 = norm(y - y_rand) / length(y);

	% Autocorrelation
	autoCorrs = CO_AutoCorr(y_rand, 1:4, 'Fourier');
	ac1 = autoCorrs(1);
	ac2 = autoCorrs(2);
	ac3 = autoCorrs(3);
	ac4 = autoCorrs(4);

	% 2-bit LZ complexity:
	% LZcomplex = EN_LZComplexity(y,3);

	% Normalized permutation entropy, PermEn(3,1) (replaced SampEn(2,0.15),
	% whose O(N^2) cost at each of the 20 randomization stages was ~75% of
	% this operation's time on long series; PermEn is O(N)):
	permEnStruct = EN_PermEn(y_rand, 3, 1);
	permen3_1 = permEnStruct.normPermEn;

	% Stationarity
	statav5 = SY_StatAv(y_rand, 'seg', 5);
	swss5_1 = SY_SlidingWindow(y_rand, 'std', 'std', 5, 1);

	out = [xcn1, xc1, d1, ac1, ac2, ac3, ac4, permen3_1, statav5, swss5_1];
end
% ------------------------------------------------------------------------------

% -------------------------------------------------------------------------------
% -------------------------------------------------------------------------------
function thehp = SUB_gethp(v)
	if v(end) > v(1)
		thehp = find(v > 0.5 * (v(end) + v(1)), 1, 'first');
	else
		thehp = find(v < 0.5 * (v(end) + v(1)), 1, 'first'); % last?
	end
end
% ------------------------------------------------------------------------------
function [c, gof] = f_fix_exp(r, dataVector, startPoint, addOffset)
	% Fits an exponential to the data vector across data points r

	% The fittype objects are built once and reused across calls (parsing the
	% model string was ~25% of each fit's cost, for ten fits per call of
	% this operation); the start point is passed to fit() instead:
	persistent fExp fExpOffset
	if isempty(fExp)
		fExp = fittype('a*exp(b*x)');
		fExpOffset = fittype('a*exp(b*x)+c');
	end
	if addOffset
		f = fExpOffset;
		f_x = @(c, x) c.a * exp(c.b * x) + c.c;
	else
		f = fExp;
		f_x = @(c, x) c.a * exp(c.b * x);
	end
	try
		[c, gof] = fit(r, dataVector, f, 'StartPoint', startPoint);
	catch
		warning('Exponential fit failed :(')
		if addOffset
			c = struct('a', NaN, 'b', NaN, 'c', NaN);
		else
			c = struct('a', NaN, 'b', NaN);
		end
		gof = struct('rsquare', NaN, 'rmse', NaN);
	end
	if doPlot;
		figure('color', 'w'); hold on;
		plot(r, dataVector, 'x-k');
		xr = linspace(min(r), max(r), 100);
		plot(xr, f_x(c, xr))
	end
end
% -------------------------------------------------------------------------------
function out = assignExpStats(out, c, gof, fieldName)
	% Assigns relevant stats from an exponential fit result, [c,gof]
	out.([fieldName, 'fexpa']) = c.a;
	out.([fieldName, 'fexpb']) = c.b;
	if (isstruct(c) && isfield(c, 'c')) || ismember('c', coeffnames(c))
		out.([fieldName, 'fexpc']) = c.c;
	end
	out.([fieldName, 'fexpr2']) = gof.rsquare;
	out.([fieldName, 'fexprmse']) = gof.rmse;
end
% -------------------------------------------------------------------------------
function out = assignExtraStats(out, dataVector, fieldName)
	% Assigns 2 extra statistics about a data vector:
	out.([fieldName, 'diff']) = abs(dataVector(end) - dataVector(1));
	out.([fieldName, 'hp']) = SUB_gethp(dataVector);
end

end
