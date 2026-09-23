function out = CO_MatrixProfile(y, m, maxN)
% CO_MatrixProfile   How well each subsequence shape recurs elsewhere in a time series.
%
% Computes the matrix profile: for every length-m window (subsequence) of the
% time series, the distance to its nearest neighbour among all other windows,
% after z-normalizing each window (so only its shape matters, not its local
% level or amplitude). Trivial matches (overlapping windows, |i-j| < m/2) are
% excluded. Distances are expressed as the equivalent nearest-neighbour Pearson
% correlation, r = 1 - d^2/(2m), which is bounded and interpretable.
%
% Features summarize the distribution of r across windows: high values mean
% shapes recur ('motifs'); a window with unusually low r is a 'discord'
% (anomaly). The corrected arc curve (Gharghabi et al., FLUSS) counts how many
% nearest-neighbour links cross each time point, relative to what a stationary
% process would give; its minimum is low when the series has a regime change,
% because windows then match within their own regime.
%
% A one-off anomaly lasting longer than about m/2 is not a discord: its own
% overlapping windows match each other (the 'twin freak' problem).
%
% Uses the STOMP recursion (O(N^2) time, O(N) memory).
%
% ---INPUTS:
% y, the input time series (z-scored in hctsa)
% m, the window length in samples; or {'ac', k} for k times the first
%       zero-crossing of the autocorrelation function (at least 10 samples).
%       Longer windows give more reliable estimates of the nearest-neighbour
%       statistics (test-retest across processes: 0.93-0.96 at k = 8 vs
%       0.73-0.84 at k = 4), at the cost of needing longer series.
% maxN, crops time series longer than this to their first maxN samples (or
%       'full' to use every sample)
%
% ---OUTPUTS:
% meanR, medianR: mean and median nearest-neighbour correlation ('matchiness')
% motifR: highest nearest-neighbour correlation (the best-repeated shape)
% discordR: lowest nearest-neighbour correlation (the most anomalous shape)
% discordGap: medianR - discordR, how anomalous the discord is relative to a
%       typical window
% propMatch90: proportion of windows with a nearest neighbour at r > 0.9
% minCAC: minimum of the corrected arc curve (low = regime change). A specialist
%       statistic: it separates regime-switching from stationary series well, but
%       across stationary series it mostly reflects estimation noise.

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
% Check inputs and set defaults:
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(m)
    m = {'ac', 8};
end
if nargin < 3 || isempty(maxN)
    maxN = 5000;
end
y = y(:);
N = length(y);
if ~(ischar(maxN) && strcmp(maxN, 'full')) && N > maxN
    y = y(1:maxN);
    N = maxN;
end

if iscell(m)
    tau = CO_FirstCrossing(y, 'ac', 0, 'discrete');
    if isnan(tau)
        out = NaN; return % data-dependent: no correlation length could be estimated
    end
    m = max(10, round(m{2}*tau));
end
numWin = N - m + 1;
exZone = ceil(m/2);
if numWin < 5*m
    out = NaN; return % data-dependent: too few windows relative to the window length
end

% ------------------------------------------------------------------------------
% Matrix profile by STOMP, in correlation units
% ------------------------------------------------------------------------------
cs = cumsum([0; y]); cs2 = cumsum([0; y.^2]);
mu = (cs(m+1:end) - cs(1:numWin))/m;
sig = sqrt(max((cs2(m+1:end) - cs2(1:numWin))/m - mu.^2, 0));
flat = sig < 1e-8*std(y);
sig(flat) = Inf; % flat windows have no shape: they match nothing (r = 0)

% Sliding dot products of the first window against all windows (by FFT):
QT = slidingDot(y(1:m), y);
QT1 = QT; % the first column, needed to restart each row of the recursion
bestR = -Inf(numWin, 1);
bestIdx = zeros(numWin, 1);
for i = 1:numWin
    if i > 1
        QT(2:end) = QT(1:end-1) - y(i-1)*y(1:numWin-1) + y(i+m-1)*y(m+1:N);
        QT(1) = QT1(i);
    end
    r = (QT - m*mu(i)*mu) ./ (m*sig(i)*sig);
    r(max(1, i-exZone+1):min(numWin, i+exZone-1)) = -Inf;
    [bestR(i), bestIdx(i)] = max(r);
end
bestR = min(bestR, 1);
bestR(flat) = NaN;
if mean(isnan(bestR)) > 0.5
    out = NaN; return % data-dependent: mostly flat windows
end

% ------------------------------------------------------------------------------
% Summaries
% ------------------------------------------------------------------------------
out.meanR = mean(bestR, 'omitnan');
out.medianR = median(bestR, 'omitnan');
out.motifR = max(bestR);
out.discordR = min(bestR);
out.discordGap = out.medianR - out.discordR;
out.propMatch90 = mean(bestR(~isnan(bestR)) > 0.9);

% Corrected arc curve: nearest-neighbour links crossing each position, relative
% to the parabola 2k(n-k)/n expected when links point to uniformly random places
ok = ~isnan(bestR);
lo = min((1:numWin)', bestIdx); hi = max((1:numWin)', bestIdx);
nc = accumarray(lo(ok), 1, [numWin+1, 1]) - accumarray(hi(ok), 1, [numWin+1, 1]);
arcs = cumsum(nc(1:numWin));
k = (1:numWin)';
ideal = 2*k.*(numWin - k)/numWin;
cac = min(arcs./ideal, 1);
edgeZone = 5*m; % the arc curve is unreliable near the edges
if numWin > 2*edgeZone
    out.minCAC = min(cac(edgeZone+1:numWin-edgeZone));
else
    out.minCAC = NaN;
end

end
% ------------------------------------------------------------------------------
function QT = slidingDot(q, t)
% Dot products of the query q with every length-m window of t
m = length(q); n = length(t);
L = 2^nextpow2(n + m);
z = real(ifft(fft(t, L) .* fft(flipud(q), L)));
QT = z(m:n);
end
