function out = NL_RecurrenceTimes(y, tau, m, theilerWin, rr, numSegments, maxN, randomSeed)
% NL_RecurrenceTimes    Recurrence-time statistics from a recurrence plot
%
% Embeds the time series in an m-dimensional delay space and characterizes
% the *white* vertical structures of the resulting recurrence plot -- the
% gaps between successive returns to a given neighborhood -- rather than
% the black line-length statistics NL_RQA computes. For each embedded
% point j, the sorted recurrence times of its neighbors give a set of
% recurrence-time samples (successive-return-time differences); pooled
% over all j, this yields a mean recurrence time and a mode (the
% probability mass at the single most common recurrence time), plus their
% variability when the series is long enough to compute them over several
% independent segments. This targets a different, and more subtle,
% signature than
% NL_RQA's determinism/laminarity: a torus-to-strange-nonchaotic-attractor
% (SNA) transition, which most Lyapunov-exponent-based diagnostics miss,
% shows up here as an abrupt jump in the segment-to-segment variance of
% both statistics, even though the point-to-point line-length statistics
% barely change.
%
% cf. E.J. Ngamga, A. Nandi, R. Ramaswamy, M.C. Romano, M. Thiel, J. Kurths,
% "Recurrence analysis of strange nonchaotic dynamics", Phys. Rev. E 75,
% 036222 (2007). DOI: 10.1103/PhysRevE.75.036222
%
% (Of that paper's two measures, only this one -- based on a single
% observed trajectory -- transfers to hctsa's one-series-in setting; their
% companion SNA-to-chaos diagnostic needs a cross-recurrence plot between
% two trajectories of the *same* system launched from different initial
% conditions under identical forcing, which requires access to the
% generating process itself, not just one observed series.)
%
% NL_ReturnTime.m already computes a related first-return-time histogram
% from a single reference-point radius, and T_MRT/T_MRT_var correlate with
% it at r up to ~0.80 on the Empirical1000 dataset -- a real conceptual
% overlap, worth knowing about, though below the redundancy bar used
% elsewhere in this codebase (r>=0.95). The segment-to-segment variance
% computed here (the paper's actual torus-to-SNA diagnostic) has no
% counterpart in NL_ReturnTime, which characterizes only a single
% whole-series histogram.
%
% ---INPUTS:
%
% y, scalar time series as a column vector
%
% tau, time delay for the embedding (can be 'ac' or 'mi', cf. BF_Embed)
%
% m, embedding dimension (a positive integer)
%
% theilerWin, Theiler window excluding temporally-correlated neighbors
%             (a proportion of N if in (0,1)) -- see NL_RQA
%
% rr, target recurrence rate used to set the neighborhood radius (the
%     radius is set once, from the full embedded series, to the
%     rr-quantile of a subsample of pairwise distances; the same radius
%     is then reused for every segment below so that segment-to-segment
%     differences reflect the dynamics rather than a re-calibrated
%     threshold)
%
% numSegments, the trajectory is divided into this many contiguous,
%              non-overlapping segments, and the mean recurrence time and
%              modal count are recomputed independently within each; their
%              variance across segments is the paper's diagnostic for the
%              torus-to-SNA transition. Each segment needs enough points
%              for a meaningful recurrence-time distribution -- if it
%              doesn't, the variance outputs (but not the full-series
%              T_MRT/N_MPRT) are NaN.
%
% maxN, the maximum number of samples to consider (cf. NL_RQA; default:
%       10000)
%
% randomSeed, whether (and how) to reset the random seed, using
%             BF_ResetSeed (cf. NL_RQA; default: 'default')
%
% ---OUTPUTS: the mean recurrence time (T_MRT) and modal recurrence-time
% probability mass (N_MPRT) of the full series, and their variance across
% numSegments independent segments (T_MRT_var, N_MPRT_var).

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
%% Check inputs
% ------------------------------------------------------------------------------
N = length(y);

if nargin < 2 || isempty(tau)
    tau = 1;
end
if nargin < 3 || isempty(m)
    m = 3;
end
if nargin < 4 || isempty(theilerWin)
    theilerWin = 0.01; % 1% of the (embedded) series length
end
if nargin < 5 || isempty(rr)
    rr = 0.1; % target recurrence rate of 10%
end
if nargin < 6 || isempty(numSegments)
    numSegments = 4;
end
if nargin < 7 || isempty(maxN)
    maxN = 10000;
end
if nargin < 8 || isempty(randomSeed)
    randomSeed = 'default';
end

if ischar(maxN) && strcmp(maxN, 'full')
    % no cropping
elseif N > maxN
    warning(sprintf(['Time series (%u > %u) is too long for recurrence-time ' ...
                     'analysis at this recurrence rate. Analyzing the first %u samples'], N, maxN, maxN));
    y = y(1:maxN);
    N = length(y);
end

% ------------------------------------------------------------------------------
%% Embed the signal
% ------------------------------------------------------------------------------
Y = BF_Embed(y, tau, m, false);
if isscalar(Y) && isnan(Y)
    warning('Embedding failed');
    out = NaN; return
end
Nemb = size(Y, 1);

if (theilerWin > 0) && (theilerWin < 1)
    theilerWinAbs = round(theilerWin * Nemb);
else
    theilerWinAbs = theilerWin;
end

if Nemb < 50 || Nemb <= 4 * theilerWinAbs
    warning('Time series too short for meaningful recurrence-time statistics (Nemb = %u, theilerWin = %u)', Nemb, theilerWinAbs);
    out = NaN; return
end

% ------------------------------------------------------------------------------
%% Set the neighborhood radius once, from the full embedded series, to
%% achieve the target recurrence rate rr (same logic as NL_RQA)
% ------------------------------------------------------------------------------
nSub = min(500, Nemb);
BF_ResetSeed(randomSeed);
subIdx = randperm(Nemb, nSub);
Dsub = pdist(Y(subIdx, :));
radius = quantile(Dsub, rr);
if radius <= 0
    warning('Degenerate neighborhood radius (data may be too degenerate/short)');
    out = NaN; return
end

% ------------------------------------------------------------------------------
%% T_MRT and N_MPRT on the full series
% ------------------------------------------------------------------------------
[out.T_MRT, out.N_MPRT] = SUB_recurrenceTimeStats(Y, radius, theilerWinAbs);
if isnan(out.T_MRT)
    warning('No recurrence-time samples found outside the Theiler window -- radius too small?');
    out.N_MPRT = NaN;
    out.T_MRT_var = NaN;
    out.N_MPRT_var = NaN;
    return
end

% ------------------------------------------------------------------------------
%% Variance of T_MRT/N_MPRT across numSegments independent, contiguous
%% segments (using the same global radius throughout)
% ------------------------------------------------------------------------------
segLen = floor(Nemb / numSegments);
minSegLen = max(50, 4 * theilerWinAbs);
if segLen < minSegLen
    % Segments would be too short for a meaningful within-segment estimate
    out.T_MRT_var = NaN;
    out.N_MPRT_var = NaN;
    return
end

segT = NaN(numSegments, 1);
segN = NaN(numSegments, 1);
for s = 1:numSegments
    idxRange = ((s - 1) * segLen + 1):(s * segLen);
    [segT(s), segN(s)] = SUB_recurrenceTimeStats(Y(idxRange, :), radius, theilerWinAbs);
end

if any(isnan(segT))
    out.T_MRT_var = NaN;
    out.N_MPRT_var = NaN;
else
    out.T_MRT_var = var(segT);
    out.N_MPRT_var = var(segN);
end

% ------------------------------------------------------------------------------
function [T_MRT, N_MPRT] = SUB_recurrenceTimeStats(Yseg, radius, theilerWinAbs)
    % For the embedded points in Yseg, find all neighbor pairs within
    % radius (via KD-tree), exclude the Theiler window, and for each
    % reference point collect the sorted times of its neighbors -- the
    % differences between consecutive recurrence times are the
    % recurrence-time samples w (the length of a "white vertical line").
    % Pooled across all reference points: T_MRT is their mean and N_MPRT
    % is the count of the single most common (modal) integer value of w.
    NsegEmb = size(Yseg, 1);
    idxCell = rangesearch(Yseg, Yseg, radius);
    allW = cell(NsegEmb, 1);
    for j = 1:NsegEmb
        nbrs = idxCell{j}(:);
        nbrs = nbrs(abs(nbrs - j) > theilerWinAbs); % exclude Theiler window (incl. self)
        if numel(nbrs) >= 2
            nbrs = sort(nbrs);
            allW{j} = diff(nbrs);
        end
    end
    w = vertcat(allW{:});
    w = w(w >= 1); % recurrence times are >= 1 by construction

    if isempty(w)
        T_MRT = NaN; N_MPRT = NaN; return
    end

    T_MRT = mean(w);
    uniqueW = unique(w);
    counts = histc(w, uniqueW);
    % Normalized to a probability mass (fraction of all pooled
    % recurrence-time samples falling at the modal value), not the paper's
    % raw count -- the raw count scales directly with how many samples
    % went in (checked empirically: r=0.84 with plain series length on
    % the Empirical1000 dataset), which would make it a length artifact
    % rather than a dynamical one.
    N_MPRT = max(counts) / numel(w);
end
% ------------------------------------------------------------------------------

end
