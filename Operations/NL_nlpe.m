function out = NL_nlpe(y, de, tau, maxN, theilerWin)
% NL_nlpe   Normalized drop-one-out constant interpolation nonlinear prediction error.
%
% Computes the nlpe for a time-delay embedded time series using Michael Small's code, nlpe
% (MS_nlpe, in Toolboxes/Michael_Small; adapted here to return the prediction residuals and
% to exclude neighbors inside a Theiler window). The series is embedded in de dimensions
% with delay tau. The one-step-ahead value following each embedded state is predicted by
% the value that followed its nearest neighbor outside the Theiler window ("drop-one-out",
% constant interpolation). The prediction errors (prediction minus data) are summarized
% using MF_ResidualAnalysis.
%
% ---INPUTS:
% y, the input time series (should be z-scored)
% de, the embedding dimension: an integer, or {'fnn',th} (or 'fnn', with th = 0.05) to
%     choose it by false nearest neighbors at threshold th, limited by the length of the
%     series analyzed: de = min(m_fnn, floor(2*log10(N))), with N the length after any
%     maxN crop. m_fnn is the smallest dimension not ruled out by NL_FNN (dimensions 1 to
%     10, a Theiler window of one autocorrelation time, escape factor 5): the first
%     dimension at which the proportion of false nearest neighbors is below th on at least
%     ceil(3/th) points, or at which too few points found a neighbor to decide, or else one
%     more than the largest dimension TISEAN computed (11 when all ten were computed); a
%     lower bound unless the proportion fell below th. The limit floor(2*log10(N)) is the
%     largest dimension N points can resolve (Eckmann & Ruelle 1992; see BF_Embed): in more
%     dimensions typical points have no neighbor within a tenth of the diameter of the set
%     (at most 4, 5, 6, 7 dimensions for N = 250, 500, 1000, 5000). Default: 3.
% tau, the time-delay (can be an integer, 'ac' for the first zero-crossing of the ACF,
%      'ac1e' for the floor of its first 1/e crossing, or 'mi' for the smaller of the
%      first minimum of the Kraskov automutual information and the 'ac1e' delay; see
%      BF_GetTau; default: 1)
% maxN, the maximum time-series length to analyze; longer series are cropped to their first
%       maxN points (default: 5000, due to memory constraints with longer time series). Set
%       to 'full' to disable, with a warning above 20000 points. The time delay, the Theiler
%       window and an 'fnn' embedding dimension (and its limit) are all computed from the
%       cropped series.
% theilerWin, Theiler window: the prediction for each point uses its nearest neighbor
%       outside this window in time ({'ac', k} for k times the first zero-crossing of the
%       autocorrelation function, or a number of samples; see BF_TheilerWindow; default:
%       {'ac', 1}). Excluding only the point itself lets a temporally adjacent point serve
%       as the neighbor, so the error measures the local smoothness of the series rather
%       than its predictability from similar past states.
%
% ---OUTPUTS: msqerr (the mean squared prediction error), and the 'full' set of
% MF_ResidualAnalysis on the prediction errors:
% meane, meanabs, stde: the mean, mean absolute value, and standard deviation of the errors
% maxonstd, the largest absolute error in units of the standard deviation of the errors
% ac1, ac2, ac3: the autocorrelation of the errors at lags 1, 2 and 3
% propbth, the proportion of the first 25 autocorrelations inside the +/-2.6/sqrt(N) band
% ftbth, the first lag at which the autocorrelation falls inside that band (26 if none)
% taurat, the ratio of the errors' first zero-crossing of the autocorrelation function to
%         that of the series
% sws, swm, the stationarity of the standard deviation and mean of the errors across 5
%         windows (SY_SlidingWindow)
% normksstat, the Kolmogorov-Smirnov statistic of the errors against a normal distribution
% popt, minsbc, the order (1 to 10) of the autoregressive model of the errors selected by
%         the Schwarz criterion, and the corresponding criterion value
% The output is NaN if tau or the embedding dimension cannot be determined, or the series
% has fewer than 20 samples.
%
% ---REFERENCES:
% M. Small, "Applied Nonlinear Time Series Analysis: Applications in Physics, Physiology,
% and Finance", World Scientific, Nonlinear Science Series A, Vol. 52 (2005).
% J.-P. Eckmann and D. Ruelle, "Fundamental limitations for estimating dimensions and
% Lyapunov exponents in dynamical systems", Physica D 56(2-3), 185-187 (1992).
% DOI: 10.1016/0167-2789(92)90023-g
%
% ---NOTES:
% Michael Small's Matlab code is available at http://small.eie.polyu.edu.hk/matlab/

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

if ~BF_iszscored(y)
	warning('The input time series should be z-scored')
end

N = length(y); % time-series length

% ------------------------------------------------------------------------------
%% Check inputs, set defaults:
% ------------------------------------------------------------------------------
% Embedding dimension:
if nargin < 2 || isempty(de)
	de = 3;
end

if nargin < 4 || isempty(maxN)
	% Set a maximum time-series length to analyze, due to memory constraints
	% with longer time series
	maxN = 5000;
end
% -------------------------------------------------------------------------------

% nlpe can cause memory pains for long time series
% Let's do this dirty cheat
if ischar(maxN) && strcmp(maxN, 'full')
	slowThreshold = 20000;
	if N > slowThreshold
		warning(['Time series (%u samples) exceeds %u with maxN=''full''; nlpe can ' ...
				 'use substantial memory for long time series'], N, slowThreshold);
	end
elseif N > maxN
	% Crop the time series to the first maxN samples:
	y = y(1:maxN);
	fprintf(1, ['Michael Small''s ''nlpe'' code is only being evaluated ' ...
				'on the first %u (/%u) time-series samples...\n'], maxN, N);
	N = maxN;
end
if N < 20 % Short time series cause problems:
	fprintf(1, 'Time series (N = %u) is too short.\n', length(y));
	out = NaN; return
end

% Time-delay, tau (from the series analyzed, i.e., after any maxN crop)
if nargin < 3 || isempty(tau)
	tau = 1;
end
if ischar(tau) && strcmp(tau, 'ac1e')
	% Adaptive delay: see BF_GetTau
	tau = BF_GetTau(y, tau);
	if isnan(tau)
		out = NaN; return
	end
end
if strcmp(tau, 'ac')
	tau = CO_FirstCrossing(y, 'ac', 0, 'discrete');
elseif strcmp(tau, 'mi')
	tau = BF_GetTau(y, 'mi');
end
if isnan(tau)
	out = NaN; return
end

% -------------------------------------------------------------------------------
% Do false nearest neighbours to compute an appropriate embedding dimension, if needed
% (escapeFactor=5 matches the convention used elsewhere for this test):
if ischar(de) || iscell(de)
	if ~iscell(de), de = {de}; end
	if ~strcmp(de{1}, 'fnn')
		error('Invalid embedding dimension specification for NL_nlpe.')
	end
	if isscalar(de)
		th = 0.05;
	else
		th = de{2};
	end
	de = NL_FNN(y, tau, 10, {'ac', 1}, 1, th, 5);
	if isnan(de)
		% NL_FNN already signals a data-dependent embedding-dimension failure
		% (e.g. series too short for the requested tau) via NaN; propagate it
		% here rather than passing a NaN de into MS_nlpe below.
		out = NaN; return
	end
	% No more dimensions than the N analyzed points can resolve (Eckmann & Ruelle 1992)
	de = min(de, max(1, floor(2*log10(N))));
end

% -------------------------------------------------------------------------------
% Run Michael Small's nonlinear prediction error code:
if nargin < 5 || isempty(theilerWin)
	theilerWin = {'ac', 1};
end
theilerWin = BF_TheilerWindow(y, theilerWin);
if isnan(theilerWin) % the autocorrelation function never crosses zero
	warning('No autocorrelation zero-crossing to set the Theiler window')
	out = NaN; return
end
res = MS_nlpe(y, de, tau, theilerWin); % residuals

% ------------------------------------------------------------------------------
%% Compute outputs
% ------------------------------------------------------------------------------
out.msqerr = mean(res.^2);
% out.rmserr = sqrt(mean(e.^2));
% out.mabserr = mean(abs(e));
% out.meanres = mean(e);

% Use MF_ResidualAnalysis on the residuals
% 1) Get statistics on residuals
residstats = MF_ResidualAnalysis(res, y, 'full');

% Convert output of residual analysis to local outputs in quick loop
fields = fieldnames(residstats);
for k = 1:length(fields)
	out.(fields{k}) = residstats.(fields{k});
end

end
