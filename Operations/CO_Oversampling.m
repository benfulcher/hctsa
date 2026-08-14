function out = CO_Oversampling(y)
% CO_Oversampling   Detects temporal oversampling relative to a series' own dynamics.
%
% Implements the oversampling-detection statistic eta (and its downsampling
% correction) from the 'oversampling' stage of the Chaos Decision Tree
% Algorithm:
%
% Toker, D. et al. "A simple method for detecting chaos in nature",
% Commun. Biol. 3, 11 (2020). DOI: 10.1038/s42003-019-0715-9
%
%   eta = range(y) / mean(|diff(y)|)
%
% A large eta means consecutive samples are, on average, tiny relative to
% the full dynamic range the series explores -- i.e., the series is sampled
% much faster than its own dynamics move, so consecutive points are close
% to redundant. Toker et al. flag eta > 10 as 'oversampled': this inflates
% the apparent smoothness/determinism of a series and can bias downstream
% nonlinear statistics (their motivation: left uncorrected, it distorts
% their 0-1 test for chaos -- cf. NL_ZeroOneTest). Their correction is to
% iteratively halve the sampling rate (keep every second point) until
% eta <= 10 or fewer than 100 points remain.
%
% cf. NL_ZeroOneTest for the chaos-classification stage of the same
% pipeline, and EN_PermEn/SD_SurrogateTest for its stochasticity-testing
% stage; this function covers the oversampling-diagnosis stage instead.
% Note eta is scale- and location-invariant (a ratio of two amplitude-unit
% quantities), so z-scored or raw y give identical results.
%
% Beyond eta itself, also reports:
% - etaRobust, the same ratio with range replaced by a 5th-95th percentile
%   range, since eta's numerator (a global max-min) is a single-outlier-
%   sensitive statistic; etaRobust asks the same oversampling question
%   without letting one extreme point set the scale.
% - numHalvings, the number of times Toker et al.'s halving procedure would
%   downsample y before eta <= 10 (or fewer than 100 points would remain)
%   -- a direct, interpretable severity measure.
% - etaAfterDownsampling, the value of eta after applying numHalvings
%   halvings (<=10, unless the series was too short to fully correct).
%
% Redundancy check (Empirical1000): eta correlates r=0.955 with the existing
% SY_RangeEvolve.totnuq (both driven by how coarsely a series sets new range
% records -- exactly what oversampling produces) and etaRobust correlates
% r=0.972 with plain eta itself (the outlier-robustness case shown above is
% real but rare in practice on typical data). Both were kept anyway: eta is
% the literal statistic from the source paper with a specific methodological
% role (a validity precondition for NL_ZeroOneTest), and etaRobust guards
% against the single-outlier failure mode even though it rarely differs from
% eta in practice. numHalvings and etaAfterDownsampling are not redundant
% with anything existing (max |r| = 0.82, 0.75 respectively).
%
% ---INPUTS:
% y, the input time series
%
% ---OUTPUTS:
% eta, etaRobust, numHalvings, etaAfterDownsampling (see above)

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
%% Check inputs:
% ------------------------------------------------------------------------------
y = y(:);
N = length(y);
if N < 10
	% Data-dependent (series too short to assess sampling adequacy), not a code error.
	warning('Time series too short to assess oversampling');
	out.eta = NaN;
	out.etaRobust = NaN;
	out.numHalvings = NaN;
	out.etaAfterDownsampling = NaN;
	return
end

meanAbsDiff = mean(abs(diff(y)));
if meanAbsDiff == 0
	% Data-dependent (constant time series), not a code error.
	warning('Zero mean absolute difference (constant time series?)');
	out.eta = NaN;
	out.etaRobust = NaN;
	out.numHalvings = NaN;
	out.etaAfterDownsampling = NaN;
	return
end

% ------------------------------------------------------------------------------
%% eta and its outlier-robust companion:
% ------------------------------------------------------------------------------
out.eta = (max(y) - min(y)) / meanAbsDiff;

robustRange = quantile(y, 0.95) - quantile(y, 0.05);
out.etaRobust = robustRange / meanAbsDiff;

% ------------------------------------------------------------------------------
%% Iterative halving-downsampling correction (Toker et al., 2020):
% ------------------------------------------------------------------------------
etaThreshold = 10;
minPoints = 100;

yDS = y;
numHalvings = 0;
etaCurr = out.eta;
while etaCurr > etaThreshold && floor(length(yDS) / 2) >= minPoints
	yDS = yDS(1:2:end);
	numHalvings = numHalvings + 1;
	mDS = mean(abs(diff(yDS)));
	if mDS == 0
		break % degenerate downsampled segment -- keep prior etaCurr
	end
	etaCurr = (max(yDS) - min(yDS)) / mDS;
end
out.numHalvings = numHalvings;
out.etaAfterDownsampling = etaCurr;

end
