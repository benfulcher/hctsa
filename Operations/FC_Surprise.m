function out = FC_Surprise(y, whatPrior, memory, numGroups, coarseGrainMethod, numIters, randomSeed)
% FC_Surprise   How surprising each next symbol is, given the recent past.
%
% Coarse-grains the time series into a sequence of symbols from a small alphabet,
% and measures how surprised a forecaster with a local memory of the past memory
% symbols would be by each new symbol. At each test point (every point with a full
% memory before it, or a fixed subsample of numIters of them for long series),
% the forecaster estimates the probability p of the symbol that actually occurred,
% using only the preceding memory symbols, and the 'information gained' (surprise)
% is -log(p), in nats.
%
% The estimated p uses Krichevsky-Trofimov-style additive smoothing,
% (numMatches + 0.5) / (nAntecedent + 0.5*alphabetSize), instead of a raw frequency
% ratio. The smoothed estimate stays strictly between 0 and 1, so the surprise is
% always finite and positive, and it falls back to the uniform guess
% 1/alphabetSize when the antecedent pattern was never observed in memory (a raw
% ratio would be undefined or 0 there, and treating "no information" as certainty
% would make a never-seen pattern the least surprising rather than the most).
%
% ---INPUTS:
% y, the input time series
%
% whatPrior, the information held in memory to predict the next symbol:
%           (i) 'dist': the distribution of symbols in the previous memory samples
%                       (default),
%           (ii) 'T1': the one-point transition probabilities in the previous
%                       memory samples (what followed the previous symbol), and
%           (iii) 'T2': the two-point transition probabilities in the previous
%                       memory samples (what followed the previous two symbols).
%
% memory, the memory length (either number of samples, or a proportion of the
%           time-series length, if between 0 and 1; default 0.2)
%
% numGroups, the number of groups to coarse-grain the time series into (default
%           3); for 'embed2quadrants' it is instead the time delay of the
%           embedding (a number of samples, 'ac1e' for the floor of the first 1/e
%           crossing of the autocorrelation function, 'mi' for the smaller of the
%           first minimum of the Kraskov automutual information and the 'ac1e'
%           delay, or 'tau' for the first zero-crossing of the autocorrelation
%           function, which is kept for backward compatibility; see BF_GetTau)
%
% coarseGrainMethod, the coarse-graining, or symbolization method (SB_CoarseGrain):
%          (i) 'quantile': an equiprobable alphabet by the value of each
%                          time-series datapoint (default),
%          (ii) 'diff': an equiprobable alphabet by the value of incremental
%                       changes in the time-series values (previously called
%                       'updown'; renamed since it is not a literal sign(diff)>0
%                       split), and
%          (iii) 'embed2quadrants': 4-letter alphabet of the quadrant each data
%                            point resides in a two-dimensional embedding space.
%
% numIters, the number of test points to use (default 500). If there are at most
%           2*numIters points with a full memory before them, all of them are used;
%           otherwise numIters of them, spread evenly by a golden-ratio sequence
%           (deterministic, and free of aliasing with periodic series).
%
% randomSeed, ignored: the test points are deterministic. Kept so that existing
%           calls still run.
%
% ---OUTPUTS:
% min, max, median, mean, sum, std: the minimum (of the nonzero values), maximum,
%       median, mean, sum and standard deviation of the surprise over the test points
% lq, uq: the lower and upper quartiles of the surprise
% propUnseen: the proportion of test points whose antecedent pattern (the current
%       symbol itself for 'dist'; the preceding 1 or 2 symbols for 'T1'/'T2') was
%       never observed anywhere in the memory window (always 0 for 'dist')
% effectSize: |mean - 1| / std, the standardized distance of the mean surprise from
%       1 nat; the length-stable form, and the one hctsa registers (NaN if the
%       surprise does not vary, std at rounding level relative to the mean)
% tstat: effectSize * sqrt(number of test points), the t-statistic of the mean
%       surprise against 1 nat (NaN under the same condition)
%
% ---NOTES:
% For 'embed2quadrants' with numGroups = 'ac1e' (what hctsa registers), the delay is
% the floor of the first 1/e crossing of the autocorrelation function (the largest
% lag at which it is still at least 1/e; see BF_GetTau), capped at floor(N/25) like
% any delay. If the autocorrelation function is undefined (a constant series) or
% never crosses 1/e, no delay exists and every output is NaN. The 1/e crossing is a
% robust measure of the correlation time; the first zero crossing (numGroups = 'tau')
% can be very long or absent (then silently replaced by the floor(N/25) cap), and is
% kept only for backward compatibility.

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
%% Check inputs and set defaults
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(whatPrior)
	whatPrior = 'dist'; % expect probabilities based on prior observed distribution
end

% memory: how far into the past to base your priors on
if nargin < 3 || isempty(memory)
	memory = 0.2; % set it as 20% of the time-series length
end
if (memory > 0) && (memory < 1) % specify memory as a proportion of the time-series length
	memory = round(memory * length(y));
end

% numGroups -- number of groups for the time-series coarse-graining/symbolization
if nargin < 4 || isempty(numGroups)
	numGroups = 3; % use three symbols to approximate the time-series values
end

% coarseGrainMethod: the coarse-graining method to use
if nargin < 5 || isempty(coarseGrainMethod)
	coarseGrainMethod = 'quantile'; % symbolize time series by their values (quantile)
end

% numIters: number of test points (all of them, unless there are more than 2*numIters)
if nargin < 6 || isempty(numIters)
	numIters = 500;
end

% ------------------------------------------------------------------------------
%% Course Grain
% ------------------------------------------------------------------------------
yth = SB_CoarseGrain(y, coarseGrainMethod, numGroups); % a coarse-grained time series using the numbers 1:numGroups

if isscalar(yth) && isnan(yth)
	% No coarse-graining exists (the embedding delay is undefined): every output is NaN
	outFields = {'min', 'max', 'median', 'mean', 'sum', 'std', 'lq', 'uq', 'propUnseen', 'effectSize', 'tstat'};
	for i = 1:length(outFields)
		out.(outFields{i}) = NaN;
	end
	return
end

% The alphabet size for the Krichevsky-Trofimov smoothing below. Usually this
% is just numGroups, but for 'embed2quadrants'/'embed2octants', numGroups is
% instead repurposed as SB_CoarseGrain's embedding time delay (numeric, or the
% string 'tau' to auto-select it) -- the alphabet size there is fixed by the
% number of quadrants/octants, not by that argument.
switch coarseGrainMethod
	case 'embed2quadrants'
		numSymbols = 4;
	case 'embed2octants'
		numSymbols = 8;
	otherwise
		numSymbols = numGroups;
end

N = length(yth); % will be the same as y, for 'quantile', and 'diff'

% Select the test points (can't test the beginning of the time series, up to memory):
numAvailable = N - memory;
if numAvailable <= 2 * numIters
	rs = (memory + 1:N)'; % every point with a full memory before it
else
	% numIters points spread over the available range by a golden-ratio sequence
	rs = unique(memory + 1 + floor(numAvailable * mod((1:numIters)' * 0.6180339887498949, 1)));
end

% -------------------------------------------------------------------------------
% Compute empirical probabilities from time series
% -------------------------------------------------------------------------------
numTest = length(rs);
store = zeros(numTest, 1); % store probabilities
nAntecedent = zeros(numTest, 1); % how many times the antecedent pattern was seen in memory
for i = 1:numTest
	switch whatPrior
		case 'dist'
			% Uses the distribution up to memory to inform the next point:
			% the "antecedent" here is trivially the whole memory window
			% (always fully observed, memory samples), so nAntecedent is
			% always memory and propUnseen will always be 0 for this prior.
			numMatches = sum(yth(rs(i) - memory:rs(i) - 1) == yth(rs(i)));
			nAntecedent(i) = memory;

		case 'T1'
			% Uses one-point correlations in memory to inform the next point

			% Estimate transition probabilities from data in memory
			% Find where in memory this has been observed before, and what
			% preceeded it:
			memoryData = yth(rs(i) - memory:rs(i) - 1);
			% Previous value observed in memory here:
			inmem = find(memoryData(1:end - 1) == yth(rs(i) - 1));
			nAntecedent(i) = length(inmem);
			if isempty(inmem)
				numMatches = 0;
			else
				numMatches = sum(memoryData(inmem + 1) == yth(rs(i)));
			end

		case 'T2'
			% Uses two-point correlations in memory to inform the next point

			memoryData = yth(rs(i) - memory:rs(i) - 1);
			% Previous value observed in memory here:
			inmem1 = find(memoryData(2:end - 1) == yth(rs(i) - 1)); % the 2:end makes the next line ok
			inmem2 = find(memoryData(inmem1) == yth(rs(i) - 2));
			nAntecedent(i) = length(inmem2);
			if isempty(inmem2)
				numMatches = 0;
			else
				% inmem2 indexes into inmem1, not directly into memoryData:
				numMatches = sum(memoryData(inmem1(inmem2) + 2) == yth(rs(i)));
			end

		otherwise
			error('Unknown method ''%s''', whatPrior);
	end
	% Krichevsky-Trofimov-style smoothed probability estimate: always in
	% (0,1), degrading to the uniform prior 1/numGroups when nAntecedent=0
	% (no information) rather than a false certainty.
	store(i) = (numMatches + 0.5) / (nAntecedent(i) + 0.5 * numSymbols);
end

% -------------------------------------------------------------------------------
% Information gained from next observation is log(1/p) = -log(p)
% -------------------------------------------------------------------------------
out.propUnseen = mean(nAntecedent == 0); % proportion of never-before-observed antecedents
store = -log(store); % transform to surprises/information gains (always finite: 0 < store < 1)
% histogram(store)

if any(store > 0)
	out.min = min(store(store > 0)); % Minimum amount of information you can gain in this way
else
	out.min = NaN;
end
out.max = max(store); % Maximum amount of information you can gain in this way
out.mean = mean(store); % mean
out.sum = sum(store); % sum (contains same information as mean since length(store) the same)
out.median = median(store); % median
out.lq = quantile(store, 0.25); % lower quartile
out.uq = quantile(store, 0.75); % upper quartile
out.std = std(store); % standard deviation

% Standardized distance of the mean information gain from 1 nat. effectSize is the
% length-stable form and is what the hctsa library registers; tstat is
% effectSize * sqrt(numTest). When the surprise does not vary (std at rounding
% level relative to the mean, e.g. for a perfectly periodic series) the ratio is
% rounding noise, so both are NaN.
if out.std < 1e-8 * out.mean
	out.effectSize = NaN;
	out.tstat = NaN;
else
	out.effectSize = abs(out.mean - 1) / out.std;
	out.tstat = out.effectSize * sqrt(numTest);
end

end
