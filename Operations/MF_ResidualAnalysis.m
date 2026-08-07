function out = MF_ResidualAnalysis(e, y, summaryLevel)
% MF_ResidualAnalysis   Canonical summary of the residuals from a model fit.
%
% This is the shared residual-summary contract for hctsa's time-series model-fitting
% and forecasting operations. Every operation that produces a residual series should
% report it through this function rather than computing its own subset, so that the
% same quantity carries the same name everywhere.
%
% Two levels are provided. 'core' is cheap (no test, no model fit) and is intended for
% operations that are themselves cheap; 'full' adds diagnostics that need extra
% machinery. The field sets were chosen by measuring 26 candidate statistics on 1,740
% residual series -- ordinary residuals from three model classes, plus GARCH
% standardized residuals -- and keeping a set in which no field is closely predictable
% from the others. 'full' recovers 94.5% of the variance of the old 26-field block.
%
% ---INPUTS:
% e, the residuals, as prediction minus data (e = yp - y), a column vector.
%
% y, (optional) the original time series the model was fitted to. Required for the
%       taurat field, which compares the residual timescale to the data timescale;
%       taurat is NaN if y is not supplied.
%
% summaryLevel, 'full' (default) or 'core'.
%
% ---OUTPUTS:
%
% core (9 fields):
%   meane      mean residual
%   meanabs    mean absolute residual (near-redundant with stde for ordinary
%              residuals, but the informative shape statistic when the residuals have
%              been standardized to unit variance, as in MF_GARCHfit)
%   stde       standard deviation of the residuals
%   maxonstd   largest absolute residual, in units of the residual standard deviation
%   ac1,ac2,ac3   residual autocorrelation at lags 1-3
%   propbth    proportion of the first 25 autocorrelations below the ~2.6/sqrt(N)
%              significance threshold
%   taurat     residual decorrelation time / data decorrelation time
%
% full (core + 6):
%   ftbth      first lag at which the autocorrelation drops below significance
%   normksstat Kolmogorov-Smirnov statistic against a normal distribution
%   sws,swm    stationarity of the residual std and mean across 5 windows
%   popt       order of an AR model fitted to the residuals (SBC-selected)
%   minsbc     the corresponding Schwarz criterion
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
%% Preliminaries
% ------------------------------------------------------------------------------
% No toolbox dependency: the spectral estimate that previously required spa() from the
% System Identification Toolbox (and, in an intermediate revision, pwelch from the Signal
% Processing Toolbox) has been removed. The four spectral fifths it produced were 83-88%
% predictable from the autocorrelation fields, which is expected -- the autocorrelation
% function and the power spectrum are Fourier duals.

if nargin < 2
	y = []; % no original series supplied; taurat will be NaN
end
if nargin < 3 || isempty(summaryLevel)
	summaryLevel = 'full';
end
if ~ismember(summaryLevel, {'core', 'full'})
	error('Unknown summary level ''%s'' (expected ''core'' or ''full'')', summaryLevel);
end

if size(e, 2) > size(e, 1)
	e = e'; % make sure residuals are a column vector
end
if all(e > 0)
	warning('Very weird that ALL model residuals are positive...')
elseif all(e < 0)
	warning('Very weird that ALL model residuals are negative...')
end
N = length(e);

% ------------------------------------------------------------------------------
%% Location, scale, and shape
% ------------------------------------------------------------------------------
out.meane = mean(e);
out.meanabs = mean(abs(e));
out.stde = std(e);
if std(e) == 0
	out.maxonstd = 0;
else
	out.maxonstd = max(abs(e)) / std(e);
end

% z-score the residuals for everything that follows (all of it is scale-free anyway):
if std(e) == 0
	eZ = zeros(N, 1);
else
	eZ = zscore(e);
end

% ------------------------------------------------------------------------------
%% Serial correlation
% ------------------------------------------------------------------------------
maxLag = 25;
autoCorrResid = CO_AutoCorr(eZ, 1:maxLag, 'Fourier');
sqrtN = sqrt(N);

out.ac1 = autoCorrResid(1);
out.ac2 = autoCorrResid(2);
out.ac3 = autoCorrResid(3);

% Proportion of the autocorrelation function within the significance band:
out.propbth = sum(abs(autoCorrResid) < 2.6 / sqrtN) / maxLag;

% Residual decorrelation time relative to that of the data. This is the most independent
% field in the block after meane -- it is not recoverable from the autocorrelations:
if isempty(y)
	out.taurat = NaN;
else
	if size(y, 2) > size(y, 1), y = y'; end
	tauY = CO_FirstCrossing(zscore(y), 'ac', 0, 'continuous');
	tauE = CO_FirstCrossing(eZ, 'ac', 0, 'continuous');
	if tauY == 0 || ~isfinite(tauY)
		out.taurat = NaN;
	else
		out.taurat = tauE / tauY;
	end
end

if strcmp(summaryLevel, 'core')
	return
end

% ------------------------------------------------------------------------------
%% (full only) Whiteness, normality, stationarity, and an AR fit to the residuals
% ------------------------------------------------------------------------------

% First lag to fall below the significance threshold (maxLag+1 if it never does):
firstBelow = find(abs(autoCorrResid) < 2.6 / sqrtN, 1, 'first');
if isempty(firstBelow)
	out.ftbth = maxLag + 1;
else
	out.ftbth = firstBelow;
end

% Normality of the residuals (the KS statistic; its p-value is a deterministic function
% of the statistic and N, and was 99.1% recoverable, so it is not reported):
[~, ~, ksstat] = kstest(eZ);
out.normksstat = ksstat;

% Are the residuals stationary? Spread of the windowed std and mean across 5 windows:
out.sws = SY_SlidingWindow(e, 'std', 'std', 5, 1);
out.swm = SY_SlidingWindow(e, 'mean', 'std', 5, 1);

% Does an AR model still find structure in the residuals?
emsg = '';
try
	[~, Aest, ~, SBC] = ARFIT_arfit(eZ, 1, 10, 'sbc', 'zero');
catch emsg
end
if ~isempty(emsg)
	warning('Error fitting AR model to residuals using ARFIT package: %s.\n', emsg.message)
	out.popt = NaN;
	out.minsbc = NaN;
else
	out.popt = length(Aest); % SBC-optimal order
	out.minsbc = min(SBC);
end

end
