function out = NL_PersistentHomology(y, tau, m, maxDim, maxN)
% NL_PersistentHomology   Topological data analysis of a time-delay embedding.
%
% Time-delay embeds the time series (Takens' theorem) and computes the
% persistent homology of the resulting point cloud under a Vietoris-Rips
% filtration, via the ripser tool (Bauer, J. Appl. Comput. Topol. 5, 391
% (2021)). A clean periodic/quasi-periodic orbit traces out a loop in
% embedding space, which shows up as one strongly persistent 1-dimensional
% (H1) topological feature; a noisy or aperiodic point cloud does not. This
% is the same underlying idea as Perea & Harer's sliding-window persistence
% approach to periodicity quantification (Found. Comput. Math. 15, 799
% (2015)), applied here directly to a plain time-delay embedding rather than
% their specific sliding-window/point-cloud construction.
%
% Persistent homology itself is standard computational topology (see, e.g.,
% Edelsbrunner & Harer, "Computational Topology: An Introduction", AMS,
% 2010); ripser computes it efficiently by never explicitly enumerating
% higher-dimensional simplices (unlike a brute-force approach, whose
% simplex-enumeration cost grows combinatorially with point-cloud size),
% which is what makes this tractable on production-length time series.
%
%---INPUTS:
% y, scalar time series as a column vector
%
% tau, time delay for the embedding (can be 'ac' or 'mi', cf. BF_Embed, or
%      'periodWelch', an option specific to this operation -- see the note in
%      the code on why a fixed tau performs much worse here than for most
%      other embedding-based operations, and the m-parameter note below for
%      why 'periodWelch' exists; default 'mi')
%
% m, embedding dimension (a positive integer; needs m >= 2 for a loop to be
%    representable at all, and m >= 3 to avoid self-intersections distorting
%    the topology of anything but the simplest orbits -- Takens' theorem
%    generically guarantees an m>=3 embedding of a 1-dimensional attractor
%    is free of such artifacts).
%
%    A companion mop is registered at a deliberately under-embedded m=2:
%    self-intersections aren't purely artifacts to be unfolded away -- a
%    curve's own shape in a *specific* low-dimensional projection can itself
%    carry information a Takens-safe embedding discards. Verified directly:
%    two-harmonic signals sin(t)+0.6*sin(2t+phi) with phi=0 vs phi=1.4 are
%    indistinguishable by any phase-blind spectral feature (|FFT| depends
%    only on harmonic amplitudes, not phi, by construction) but separate
%    cleanly at m=2 -- almost entirely as a difference in the single
%    dominant loop's persistence, not its count, despite the naive
%    expectation of a literal figure-eight (two independent 1-cycles); see
%    the m=3 default, where the gap collapses substantially, consistent
%    with the effect being tied to the specific 2D projection rather than
%    the intrinsic topology of the (here, truly 1-dimensional) attractor.
%
%    This m=2 companion mop uses tau='periodWelch' (see the tau-note in the
%    code), not 'mi'. 'mi' was tried first and looked promising (maxPersistenceH1
%    0.35 vs 0.63, ~30 std devs apart across 30 reps) but turned out to be
%    largely lucky, not principled: a tau sweep over the same construction
%    showed the discrimination effect's sign/magnitude oscillates with the
%    embedding delay (ratio ranged 0.12x-1.81x across integer tau), and 'mi'
%    happened to land on the single best point in that sweep at only one of
%    three sampling rates tested -- at the other two, it agreed with the
%    (much worse-performing) 'ac' heuristic. A harmonic-ratio variant of the
%    same construction (1:3 instead of 1:2) additionally showed the *sign*
%    of the effect isn't stable across related constructions either.
%    'periodWelch' doesn't fix the sign-instability (nothing about tau
%    selection can), but it does fix the sampling-rate sensitivity: it
%    targets a fixed *physical* delay (not an absolute index count) by
%    normalizing against the series' own estimated dominant period, so the
%    same relative embedding delay is used regardless of how finely the
%    series happens to be sampled.
%
% maxDim, maximum homology dimension to compute (0 and 1 are the only
%         dimensions used by the outputs below; higher dimensions are far
%         more expensive and not currently exposed as separate fields)
%
% maxN, the maximum number of embedded points to feed to ripser. Vietoris-
%       Rips persistent homology cost grows steeply with point-cloud size, so
%       longer embeddings are reduced to maxN points via an evenly-spaced
%       (not random, not a truncation to an early window) subsample of the
%       embedded point cloud -- this keeps coverage of the full reconstructed
%       attractor rather than just an early time window, unlike e.g. NL_RQA's
%       maxN, which crops the raw series instead (RQA's per-pair cost is what
%       forces that choice there; here the embedding itself is cheap and only
%       the point count handed to ripser needs capping).
%
%---OUTPUTS: summary statistics of the H1 (loop) persistence diagram --
% maximum persistence (strength of the single most persistent loop), total
% persistence, and the Shannon entropy of the persistence distribution (one
% dominant loop vs. many comparable ones) -- plus total H0 persistence,
% capturing clustering structure of the point cloud at small scales. All
% persistence values are normalized by the Rips filtration threshold used,
% so they are comparable in scale across different input series.
%
% (An earlier candidate field, a count of "significant" H1 intervals above a
% fixed normalized-persistence threshold, was dropped after checking against
% real data (Bonn EEG, Empirical1000): it was 0 for 149/150 Bonn EEG series
% -- real (noisy, imperfectly periodic) data essentially never clears a
% fixed significance bar tuned by eye against a clean sine wave, so the
% field carried almost no information. maxPersistenceH1 above captures the
% same underlying signal as a continuous quantity instead.)

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
if nargin < 2 || isempty(tau)
    % 'mi' (first minimum of mutual information), not a fixed tau=1: loop
    % shape -- and hence H1 persistence relative to point-cloud diameter --
    % is highly sensitive to this choice. A fixed tau=1 embedding of a
    % smooth periodic series is nearly collinear locally (adjacent samples
    % are highly correlated), producing a thin, elongated loop whose
    % persistence is small relative to the cloud's overall diameter;
    % validated directly (clean sine wave, tau=1): loop persistence only
    % ~14% of the diameter, indistinguishable from noise. The same series
    % with tau='mi' gives ~55% -- an order of magnitude cleaner separation
    % from noise (~4%) and noisy/chaotic series (~1-6%).
    tau = 'mi';
end
if nargin < 3 || isempty(m)
    m = 3;
end
if nargin < 4 || isempty(maxDim)
    maxDim = 1;
end
if nargin < 5 || isempty(maxN)
    % Cost grows steeply with point-cloud size (measured on this machine, a
    % noisy-periodic series capped to N points, dim<=1): ~0.14s at 500,
    % ~0.56s at 1000, ~3.3s at 2000, ~15s at 4000. 1000 is a reasonable
    % default balance of signal quality vs. per-series cost at hctsa's scale.
    maxN = 1000;
end

% ------------------------------------------------------------------------------
%% Embed the signal
% ------------------------------------------------------------------------------
if ischar(tau) && strcmpi(tau, 'periodWelch')
    % BF_Embed doesn't know this option -- resolve it to a plain number (or
    % fall back to the string 'mi', which BF_Embed does know how to resolve,
    % including its own NaN-safety for too-short series) before embedding.
    tau = SUB_periodNormalizedTau(y);
end
X = BF_Embed(y, tau, m, false);
if isscalar(X) && isnan(X) % embedding failed
    out = NaN; return
end
Nemb = size(X, 1);

if Nemb < 20
    % Too few embedded points for a meaningful Rips filtration
    out = NaN; return
end

if Nemb > maxN
    % Evenly-spaced subsample over the full embedded point cloud (see header
    % comment for why this differs from a simple crop-to-first-maxN):
    idx = round(linspace(1, Nemb, maxN));
    X = X(idx, :);
    Nemb = size(X, 1);
end

% ------------------------------------------------------------------------------
%% Estimate a Rips filtration threshold (point-cloud diameter) from a cheap
%% deterministic subsample, so ripser isn't asked to build the filtration
%% all the way up to a (possibly much larger, outlier-driven) true diameter
% ------------------------------------------------------------------------------
nSub = min(300, Nemb);
subIdx = round(linspace(1, Nemb, nSub));
Dsub = pdist(X(subIdx, :));
threshold = max(Dsub);
if ~(threshold > 0)
    % Degenerate point cloud (e.g., constant series)
    out = NaN; return
end

% ------------------------------------------------------------------------------
%% Run ripser
% ------------------------------------------------------------------------------
filePath = BF_WriteTempFile(X);
ripser_command = sprintf('ripser --format point-cloud --dim %u --threshold %.8g %s', ...
                            maxDim, threshold, filePath);
[~, res] = BF_RipserSystem(ripser_command);

H1 = SUB_parseDim(res, 1);
H0 = SUB_parseDim(res, 0);

% ------------------------------------------------------------------------------
%% Summarize the H1 (loop) persistence diagram
% ------------------------------------------------------------------------------
if ~isempty(H1)
    pers1 = H1(:, 2) - H1(:, 1);
    pers1 = pers1(isfinite(pers1)) / threshold; % normalize; drop essential class(es)
else
    pers1 = [];
end

if isempty(pers1)
    out.maxPersistenceH1 = 0;
    out.totalPersistenceH1 = 0;
    out.persistenceEntropyH1 = 0;
else
    out.maxPersistenceH1 = max(pers1);
    out.totalPersistenceH1 = sum(pers1);
    % Dense Rips complexes routinely produce exact-zero-persistence H1 pairs
    % (simplices born and paired off at the same filtration value); a p==0
    % term must contribute 0 to the entropy sum by convention, but MATLAB's
    % 0*log(0) evaluates to NaN (an indeterminate form), silently corrupting
    % the whole sum. Drop zero-persistence entries before computing p so
    % this can't happen (found via a clean periodic signal that produced
    % exactly this case).
    pers1pos = pers1(pers1 > 0);
    p = pers1pos / sum(pers1pos);
    out.persistenceEntropyH1 = -sum(p .* log(p));
end

% ------------------------------------------------------------------------------
%% H0 summary (clustering structure of the point cloud)
% ------------------------------------------------------------------------------
if ~isempty(H0)
    pers0 = H0(:, 2) - H0(:, 1);
    pers0 = pers0(isfinite(pers0)) / threshold;
else
    pers0 = [];
end
out.totalPersistenceH0 = sum(pers0); % 0 if empty

% ------------------------------------------------------------------------------
function tau = SUB_periodNormalizedTau(y)
    % Targets a fixed *physical* embedding delay (5, in units of the
    % series' own estimated dominant period / 2*pi -- see the m-parameter
    % docstring note for why this specific value) rather than an absolute
    % index count, so the resulting index-tau automatically scales with
    % however finely the series happens to be sampled.
    %
    % Two failure modes ruled out during development, in order:
    % (1) a plain single-FFT argmax period estimate locks onto low-frequency
    %     trend/drift on real (broadband) data -- median "period" of 100+
    %     samples on Bonn EEG/Empirical1000 test series, with some series
    %     maxing out at half the series length (i.e., just the DC-adjacent
    %     bin). Real time series rarely have one clean dominant sinusoid the
    %     way a synthetic test construction does.
    % (2) findpeaks on a single (non-averaged) periodogram fixes (1) but has
    %     essentially zero specificity: it "finds" a locally-prominent peak
    %     on 100% of pure-white-noise realizations tested, at every
    %     prominence threshold tried, because periodogram bin-to-bin
    %     variance alone produces apparent local peaks by chance.
    % Welch's method (averaged periodogram over overlapping segments)
    % suppresses that noise-driven variance directly. Validated on synthetic
    % data: at MinPeakProminence=2 (log scale, ~7.4x), 0% false positives on
    % near-unit-root AR(1) (trend, no periodicity), 10% on white noise,
    % 100% true-positive rate with accurate period recovery on a sine wave
    % even under heavy added noise. On real data (Bonn EEG, Empirical1000),
    % a significant peak was found for ~80% of series, with resulting tau
    % values comparable in scale to 'mi''s own (median ~8-12 vs 'mi''s
    % median ~4-10) -- unlike the two ruled-out methods above, whose
    % resulting tau values were wildly unstable (medians 20-90+, maxima in
    % the hundreds to thousands).
    targetPhysDelay = 5;
    minProm = 2.0;

    if numel(y) < 16
        % pwelch errors outright (rather than degrading gracefully) below
        % ~8 samples with default windowing; treat any series this short as
        % too short for a meaningful period estimate regardless.
        tau = 'mi'; return
    end
    yDetrended = y - mean(y);
    [Pxx, F] = pwelch(yDetrended, [], [], [], 1);
    Pxx = Pxx(2:end); F = F(2:end); % exclude the DC bin
    if numel(Pxx) < 3
        tau = 'mi'; return
    end
    logP = log(Pxx + eps(max(Pxx)));
    [~, locs, ~, proms] = findpeaks(logP, 'MinPeakProminence', minProm);
    if isempty(locs)
        % No sufficiently prominent period found -- period-normalization
        % isn't meaningful without one; fall back to the general-purpose
        % heuristic (about 20% of real series in validation).
        tau = 'mi'; return
    end
    [~, best] = max(proms); % most prominent peak, not necessarily tallest
    period = 1 / F(locs(best));
    tau = max(1, round(period * targetPhysDelay / (2*pi)));
end
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
function ivals = SUB_parseDim(res, d)
    % Extracts [birth, death] pairs for dimension d from ripser's stdout;
    % an unbounded ("essential") interval, printed with a blank death field
    % (' [b, )'), is returned with death = Inf.
    pat = sprintf('persistence intervals in dim %d:\\n((?:.*\\n)*?)(?=persistence intervals in dim|$)', d);
    tok = regexp(res, pat, 'tokens', 'once');
    if isempty(tok)
        ivals = zeros(0, 2); return
    end
    pairs = regexp(tok{1}, '\[([^,]*),([^\)]*)\)', 'tokens');
    ivals = zeros(numel(pairs), 2);
    for i = 1:numel(pairs)
        birthStr = strtrim(pairs{i}{1});
        deathStr = strtrim(pairs{i}{2});
        ivals(i, 1) = str2double(birthStr);
        if isempty(deathStr)
            ivals(i, 2) = Inf;
        else
            ivals(i, 2) = str2double(deathStr);
        end
    end
end
% ------------------------------------------------------------------------------

end
