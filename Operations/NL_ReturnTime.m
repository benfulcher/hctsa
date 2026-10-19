function out = NL_ReturnTime(y, NNR, numLags, past, Nref, embedParams)
% NL_ReturnTime   Analysis of the histogram of return times.
%
% Return times are the time taken for the time series to return to a similar location in
% phase space for a given reference point. Strong peaks in the histogram are indicative of
% periodicities in the data.
%
% For each reference point in the embedding space, its NNR nearest neighbors are found
% (excluding a Theiler window of "past" samples either side), and the time offset, T, of
% each neighbor from the reference point is recorded. The histogram of these offsets over
% the numLags lags beyond the Theiler window, T = past+1, ..., past+numLags (the "return-time
% profile"), is analyzed. This follows TSTOOL's 'return_time', with one change: each
% lag's count is divided by its expected count if neighbors
% were placed at random among the valid (Theiler-excluded) candidates, rather than TSTOOL's
% 2*NNR*(N - T), so that the histogram is ~1 at every lag for an uncorrelated process at
% any series length (TSTOOL's normalization scaled as 1/N). Values above 1 mark lags at
% which the trajectory preferentially returns to its neighborhood. The profile is closely
% related to the tau-recurrence rate of recurrence quantification analysis, with
% neighborhoods holding a fixed proportion of points rather than having a fixed radius.
%
% The neighborhood of a reference point is the ball out to its NNR-th nearest neighbor,
% boundary included, with distances that differ by at most 1e-12 of the range of y treated
% as equal: so the neighborhood, and the profile, do not depend on rounding or on the order
% in which neighbors are found. With exactly repeated delay vectors (quantized data,
% periodic orbits) a neighborhood can hold more than NNR points, and the expected count
% uses the number it holds.
%
% (For the distribution of *first* return times to a neighborhood, see NL_RecurrenceTimes.)
%
% ---INPUTS:
% y, scalar time series as a column vector
% NNR, number of nearest neighbors (or, if in (0,1), a proportion of the number of embedded
%      points, keeping neighborhoods the same size in probability as the series length
%      changes; default: 0.01)
% numLags, the number of lags beyond the Theiler window to analyze, in samples (default: 100)
% past, Theiler window, excluding neighbors that are close only because they are close in
%       time: {'ac', k} for k times the first zero-crossing of the autocorrelation
%       function, or a number of samples (see BF_TheilerWindow; default: {'ac', 1})
% Nref, number of reference points, spaced evenly through the series (-1 uses all points;
%       default: -1). A fixed number keeps the number of neighbors counted at each lag, and
%       so the sampling noise of the histogram, independent of the series length (neighbors
%       are still sought among all points).
% embedParams, the embedding, as {tau, m}, to feed into BF_Embed (default: {'ac', 'fnn'})
%
% ---OUTPUTS: measures of the return-time profile (the neighbor count at each lag relative to
% chance), and of the histogram of its values:
% max, std, iqr, the maximum, standard deviation and interquartile range of the profile
% pzeros, the proportion of lags with no neighbors
% pg05, the proportion of lags at which the profile exceeds half its maximum
% meanpeaksep, maxpeaksep, minpeaksep, rangepeaksep, stdpeaksep, statistics of the intervals
%      between successive crossings of half the maximum (which are alternately the widths of
%      peaks and the gaps between them, those shorter than 0.4 of the longest interval being
%      dropped), as a proportion of the number of lags (stdpeaksep is divided by the square
%      root of the number of lags instead; all are NaN with fewer than 3 crossings)
% fracsdrt, fracmeanrt, the share of the first half of the lags in the standard deviation
%      (fracsdrt) or mean (fracmeanrt) of the profile: a/(a + b), with a and b the statistic
%      over the first and second halves. This lies in [0,1]: 0.5 when the halves match and 1
%      when all the variation (or mass) is at short lags; NaN if both halves are constant
%      (or empty)
% hhist, the entropy of the profile as a distribution over lags
% hcgdist, rangecgdist, pzeroscgdist, the entropy, range and proportion of zeros of the
%      profile after summing it into 20 equal bins of lags (as a distribution over bins)
% maxhisthist, phisthistmin, hhisthist, the maximum, the first (lowest-value) bin
%      probability, and the entropy of the histogram of profile values (ceil(sqrt(n)) equal-width bins, BF_HistEdges)
%
% ---REFERENCES:
% Marwan et al., Phys. Rep. 438, 237 (2007) (recurrence quantification analysis).

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
%% Check Inputs
% ------------------------------------------------------------------------------
% Number of nearest neighbours, NNR (a proportion is resolved after embedding)
if nargin < 2 || isempty(NNR)
	NNR = 0.01;
end

% Number of lags to analyze, numLags
if nargin < 3 || isempty(numLags)
	numLags = 100;
end
if numLags < 2
	error('numLags (%g) must be at least 2', numLags);
end

% Theiler window, past
if nargin < 4 || isempty(past)
	past = {'ac', 1};
end
past = BF_TheilerWindow(y, past);
if isnan(past) % the autocorrelation function never crosses zero
	warning('No autocorrelation zero-crossing to set the Theiler window')
	out = NaN; return
end
maxT = past + numLags; % maximum return time (lag) to consider

% Number of reference points
if nargin < 5 || isempty(Nref)
	Nref = -1; % use all available points
end

% embed parameters
if nargin < 6 || isempty(embedParams)
	embedParams = {'ac', 'fnn'};
	fprintf(1, 'Using default embedding using autocorrelation and false nearest neighbors\n');
end

doPlot = false; % plot outputs to figures

% ------------------------------------------------------------------------------
%% Embed the signal (native MATLAB matrix embedding, not a TSTOOL/TISEAN call)
% ------------------------------------------------------------------------------
Y = BF_Embed(y, embedParams{1}, embedParams{2}, false);
if isscalar(Y) && isnan(Y) % embedding failed
	warning('Embedding failed');
	out = NaN; return
end
N_embed = size(Y, 1);
if (NNR > 0) && (NNR < 1) % a proportion of the number of embedded points
	NNR = max(1, round(NNR * N_embed));
end
if N_embed < 2 * maxT || N_embed <= NNR + 2 * past + 1
	% Need every lag in the histogram to be sampled by at least half the points
	warning('Time series too short to do a return-time analysis with these parameters')
	out = NaN; return
end

% ------------------------------------------------------------------------------
%% Neighborhood radius of each reference point: distance to its NNR-th nearest
%% neighbor outside the Theiler window
% ------------------------------------------------------------------------------
if Nref == -1 || Nref >= N_embed
	refIdx = (1:N_embed)';
else
	refIdx = unique(round(linspace(1, N_embed, Nref)))';
end
% At most 2*past + 1 points (the reference point itself included) fall within
% the Theiler window, so NNR + 2*past + 1 neighbors always hold NNR valid ones
K = min(N_embed, NNR + 2 * past + 1);
tol = 1e-12 * (max(y) - min(y)); % distances closer than this are equal
r2 = NaN(N_embed, 1); % squared radius (NaN for non-reference points)
chunkSize = max(1, floor(2e6 / K)); % bound the memory of the neighbor lists
for c = 1:chunkSize:length(refIdx)
	theRefs = refIdx(c:min(c + chunkSize - 1, length(refIdx)));
	[idx, dist] = knnsearch(Y, Y(theRefs, :), 'K', K);
	isValid = abs(idx - theRefs) > past;
	[~, whichCol] = max(cumsum(isValid, 2) >= NNR, [], 2);
	% (the boundary is included, with the tolerance, so that every point at the
	% NNR-th neighbor's distance is in, whatever the arithmetic path to it)
	r2(theRefs) = (dist(sub2ind(size(dist), (1:length(theRefs))', whichCol)) + tol).^2;
end

% The number of valid (Theiler-excluded) points in each neighborhood: NNR, or
% more where points tie with the NNR-th neighbor
nNear = zeros(N_embed, 1);
chunkSize = max(1, floor(2e6 / N_embed));
for c = 1:chunkSize:length(refIdx)
	theRefs = refIdx(c:min(c + chunkSize - 1, length(refIdx)))';
	D2 = zeros(N_embed, length(theRefs));
	for j = 1:size(Y, 2)
		D2 = D2 + (Y(:, j) - Y(theRefs, j)').^2;
	end
	inBall = (D2 <= r2(theRefs)') & (abs((1:N_embed)' - theRefs) > past);
	nNear(theRefs) = sum(inBall, 1);
end

% ------------------------------------------------------------------------------
%% Count neighbors at each lag, relative to the count expected by chance
% ------------------------------------------------------------------------------
lags = (past + 1:maxT)';
numLags = length(lags);
counts = zeros(numLags, 1);
for k = 1:numLags
	T = lags(k);
	fwd = refIdx(refIdx + T <= N_embed); % references with a partner T ahead
	bwd = refIdx(refIdx - T >= 1); % references with a partner T behind
	counts(k) = sum(sum((Y(fwd + T, :) - Y(fwd, :)).^2, 2) <= r2(fwd)) ...
			+ sum(sum((Y(bwd - T, :) - Y(bwd, :)).^2, 2) <= r2(bwd));
end
% By chance, a given valid candidate is one of the n_i points in reference i's
% neighborhood with probability n_i/V_i, where V_i is the number of points outside
% i's Theiler window
i = (1:N_embed)';
V = N_embed - (min(i - 1, past) + min(N_embed - i, past) + 1);
w = zeros(N_embed, 1);
w(refIdx) = nNear(refIdx) ./ V(refIdx);
cw = cumsum(w);
expected = cw(N_embed - lags) + (cw(end) - cw(lags)); % forward + backward partners
Trett = counts ./ expected;

if doPlot
	figure('color', 'w'); %#ok<UNRCH>
	plot(lags, Trett, 'k')
	xlabel('Lag, T'); ylabel('Neighbors relative to chance')
end

% ------------------------------------------------------------------------------
%% Quantify structure in output
% ------------------------------------------------------------------------------
NN = numLags;
out.max = max(Trett);
out.std = std(Trett);
out.pzeros = sum(Trett == 0) / NN;
out.pg05 = sum(Trett > max(Trett) * 0.5) / NN;
out.iqr = iqr(Trett);

% recurrent peaks:
icross05 = find((Trett(1:end - 1) - 0.5 * max(Trett)) .* (Trett(2:end) - 0.5 * max(Trett)) < 0);
if ~isempty(icross05) && length(icross05) > 2
	difficross05 = diff(icross05);
	difficross05 = difficross05(difficross05 > 0.4 * max(difficross05)); % remove small entries, crossing peaks

	out.meanpeaksep = mean(difficross05) / NN;
	out.maxpeaksep = max(difficross05) / NN;
	out.minpeaksep = min(difficross05) / NN;
	out.rangepeaksep = range(difficross05) / NN;
	out.stdpeaksep = std(difficross05) / sqrt(NN);
else
	out.meanpeaksep = NaN;
	out.maxpeaksep = NaN;
	out.minpeaksep = NaN;
	out.rangepeaksep = NaN;
	out.stdpeaksep = NaN;
end

% short lags compared to long lags, as the first half's share of each statistic
% (bounded, and defined unless both halves are constant, or empty):
tolTrett = 1e-12 * max(Trett);
sdShort = std(Trett(1:floor(end / 2)));
sdLong = std(Trett(floor(end / 2) + 1:end));
meanShort = mean(Trett(1:floor(end / 2)));
meanLong = mean(Trett(floor(end / 2) + 1:end));
if sdShort + sdLong > tolTrett
	out.fracsdrt = sdShort / (sdShort + sdLong);
else
	out.fracsdrt = NaN;
end
if meanShort + meanLong > tolTrett
	out.fracmeanrt = meanShort / (meanShort + meanLong);
else
	out.fracmeanrt = NaN;
end

% entropy of the histogram, as a distribution over lags:
pTrett = Trett / sum(Trett);
out.hhist = -sum(pTrett(pTrett > 0) .* log(pTrett(pTrett > 0)));

% ------------------------------------------------------------------------------
%% Coarse-grain to 20 bins of lags
% ------------------------------------------------------------------------------
numBins = 20;
cglav = zeros(numBins, 1);
inds = round(linspace(0, NN, numBins + 1));
for i = 1:numBins
	cglav(i) = sum(pTrett(inds(i) + 1:inds(i + 1)));
end
if doPlot
	figure('color', 'w'); %#ok<UNRCH>
	box('on');
	plot(cglav, 'k')
end
out.hcgdist = -sum(cglav(cglav > 0) .* log(cglav(cglav > 0)));
out.rangecgdist = range(cglav);
out.pzeroscgdist = sum(cglav == 0) / numBins;

% ------------------------------------------------------------------------------
%% Get distribution of distribution of return times
% ------------------------------------------------------------------------------
binEdges = BF_HistEdges(Trett, 'sqrt'); % ceil(sqrt(n)) equal-width bins spanning the values
nhist = histcounts(Trett, binEdges, 'Normalization', 'probability');
if doPlot
	binCenters = mean([binEdges(1:end - 1); binEdges(2:end)]); %#ok<UNRCH>
	figure('color', 'w');
	plot(binCenters, nhist, 'o-k')
end
out.maxhisthist = max(nhist);
out.phisthistmin = nhist(1); % probability in the first (lowest-value) bin
out.hhisthist = -sum(nhist(nhist > 0) .* log(nhist(nhist > 0)));

end
