function out = NL_ReturnTime(y, NNR, maxT, past, Nref, embedParams)
% NL_ReturnTime    Analysis of the histogram of return times.
%
% Return times are the time taken for the time series to return to a similar
% location in phase space for a given reference point.
%
% Strong peaks in the histogram are indicative of periodicities in the data.
%
% ---INPUTS:
%
% y, scalar time series as a column vector
% NNR, number of nearest neighbours
% maxT, maximum return time to consider
% past, Theiler window
% Nref, number of reference indicies
% embedParams, to feed into BF_Embed
%
% ---OUTPUTS: include basic measures from the histogram, including the occurrence of
% peaks, spread, proportion of zeros, and the distributional entropy.
%
% Computed natively in MATLAB (this operation previously used TSTOOL's
% 'return_time'; TISEAN has no direct equivalent -- its own recurrence
% tool, 'recurr', defines neighborhoods by a fixed epsilon radius rather
% than a nearest-neighbor count, a parameter-type mismatch with NNR below).
% For each of Nref reference points, the neighborhood radius is the
% distance to its NNR-th nearest neighbor (excluding a Theiler window of
% "past" samples, via a KD-tree, same approach as NL_localdensity.m);
% starting just after that Theiler window, the series is scanned forward
% for the first return within that radius, up to maxT samples ahead. A
% return time of 0 is a sentinel for "no return found within maxT" (a
% return time can never be legitimately 0, since the Theiler window
% already excludes any offset from 0 up to "past").

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
N = length(y); % length of the input time series

% Number of nearest neighbours, NNR
if nargin < 2 || isempty(NNR)
	NNR = 5;
end
if (NNR > 0) && (NNR < 1) % specify a proportion of time series length
	NNR = floor(NNR * N); if NNR == 0, NNR = 1; end
end

% Maximum return time, maxT
if nargin < 3 || isempty(maxT)
	maxT = 0.1;
end
if (maxT > 0) && (maxT <= 1) % specify a proportion
	maxT = floor(N * maxT);
	if maxT == 0, maxT = 1; end
end

% Theiler window, past
if nargin < 4 || isempty(past)
	past = 10;
end
if (past > 0) && (past < 1) % specify a proportion
	past = floor(N * past);
	if past == 0, past = 1; end % round up from 0
end

% Number of reference points
if nargin < 5 || isempty(Nref)
	Nref = -1; % use all available points
end

% embed parameters
if nargin < 6 || isempty(embedParams)
	embedParams = {'ac', 'fnn'};
	fprintf(1, 'Using default embedding using autocorrelation and cao\n');
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
[N_embed, m] = size(Y);
if N_embed < 10
	% Set heuristic minimum (10) on the number of points needed to perform a meaningful analysis
	warning('Time series not long enough for return time analysis')
	out = NaN; return
end
if N_embed <= NNR + 2 * past
	warning('Time series too short to do a return-time analysis with these parameters')
	out = NaN; return
end

% ------------------------------------------------------------------------------
%% Resolve the reference points
% ------------------------------------------------------------------------------
if Nref == -1 || Nref >= N_embed
	refIdx = 1:N_embed;
else
	refIdx = randperm(N_embed, Nref); % a random subset of reference points
end
NN = length(refIdx);

% ------------------------------------------------------------------------------
%% For each reference point, find its NNR-th-nearest-neighbor radius (KD-tree,
%% Theiler-window-aware, same approach as NL_localdensity.m), then scan
%% forward for the first return within that radius
% ------------------------------------------------------------------------------
kFetch = min(N_embed - 1, NNR + 2 * past + 5);
[idx, dist] = knnsearch(Y, Y(refIdx, :), 'K', kFetch + 1);

Trett = zeros(NN, 1);
for ii = 1:NN
	i = refIdx(ii);

	validDists = dist(ii, abs(idx(ii, :) - i) > past);
	if length(validDists) < NNR
		allDists = sqrt(sum((Y - Y(i, :)).^2, 2));
		allDists(abs((1:N_embed)' - i) <= past) = Inf;
		validDists = sort(allDists);
	end
	r_i = validDists(NNR);

	winStart = i + past + 1;
	winEnd = min(i + maxT, N_embed);
	if winStart > winEnd
		Trett(ii) = 0; % sentinel: no return found within maxT
		continue
	end
	candidateIdx = winStart:winEnd;
	sqDists = sum((Y(candidateIdx, :) - Y(i, :)).^2, 2);
	firstReturn = find(sqDists <= r_i^2, 1, 'first');
	if isempty(firstReturn)
		Trett(ii) = 0; % sentinel: no return found within maxT
	else
		Trett(ii) = candidateIdx(firstReturn) - i;
	end
end

% ------------------------------------------------------------------------------
%% Quantify structure in output
% ------------------------------------------------------------------------------
out.max = max(Trett);
out.std = std(Trett(Trett > 0)); % exclude the "no return found" sentinel
out.pzeros = sum(Trett == 0) / NN;
out.pg05 = sum(Trett > max(Trett) * 0.5) / NN;
out.iqr = iqr(Trett(Trett > 0)); % exclude the "no return found" sentinel

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

% exclude the "no return found" sentinel from each half before comparing:
TrettFirstHalf = Trett(1:floor(end / 2));
TrettFirstHalf = TrettFirstHalf(TrettFirstHalf > 0);
TrettSecondHalf = Trett(floor(end / 2) + 1:end);
TrettSecondHalf = TrettSecondHalf(TrettSecondHalf > 0);
out.statrtys = std(TrettFirstHalf) / std(TrettSecondHalf);
out.statrtym = mean(TrettFirstHalf) / mean(TrettSecondHalf);

out.hhist = -sum(Trett(Trett > 0) .* log(Trett(Trett > 0)));

% ------------------------------------------------------------------------------
%% Coarse-grain to 20 bins
% ------------------------------------------------------------------------------
numBins = 20;
cglav = zeros(numBins, 1);
inds = round(linspace(0, NN, numBins + 1));
for i = 1:numBins
	cglav(i) = sum(Trett(inds(i) + 1:inds(i + 1)));
end
if doPlot
	figure('color', 'w');
	box('on');
	plot(cglav, 'k')
end
out.hcgdist = -sum(cglav(cglav > 0) .* log(cglav(cglav > 0)));
out.rangecgdist = range(cglav);
out.pzeroscgdist = sum(cglav == 0) / numBins;

% ------------------------------------------------------------------------------
%% Get distribution of distribution of return times
% ------------------------------------------------------------------------------
[nhist, binEdges] = histcounts(Trett, 'BinMethod', 'sqrt', 'Normalization', 'probability');
if doPlot
	binCenters = mean([binEdges(1:end - 1); binEdges(2:end)]);
	figure('color', 'w');
	plot(binCenters, nhist, 'o-k')
end
out.maxhisthist = max(nhist);
out.phisthistmin = nhist(1); % probability in the first (smallest-return-time) bin
out.hhisthist = -sum(nhist(nhist > 0) .* log(nhist(nhist > 0)));

end
