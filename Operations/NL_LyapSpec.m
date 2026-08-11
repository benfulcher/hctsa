function out = NL_LyapSpec(y, tauMethod, m, kNN, maxN)
% NL_LyapSpec   Spectrum of Lyapunov exponents of a time series.
%
% Estimates the full spectrum of m Lyapunov exponents (not just the
% largest -- cf. NL_LargestLyap, which wraps TISEAN's 'lyap_r', a
% single-nearest-neighbor-tracking method that only ever returns one
% number) using TISEAN's 'lyap_spec', an implementation of the
% Sano-Sawada method: a time-delay embedding of y is built, and local
% Jacobians are estimated from each point's neighbors and propagated
% forward with periodic QR reorthogonalization to give m exponents,
% ordered from largest to smallest.
%
% This is a qualitatively different kind of information to
% NL_LargestLyap: rather than a single "is it chaotic and how fast does
% it diverge" number, the full spectrum gives the shape of local
% expansion/contraction across all m embedded directions, from which
% overall dissipation rate (sum of all exponents) and an estimate of the
% attractor's fractal dimension (Kaplan-Yorke conjecture) can both be
% derived, cf. ---OUTPUTS below.
%
% cf. M. Sano & Y. Sawada, "Measurement of the Lyapunov spectrum from a
% chaotic time series", Phys. Rev. Lett. 55(10), 1082 (1985).
%
% J. Kaplan & J. Yorke, "Chaotic behavior of multidimensional difference
% equations", in Functional Differential Equations and Approximation of
% Fixed Points, Lecture Notes in Mathematics 730, 204-227 (1979) -- for
% the Kaplan-Yorke dimension conjecture used in KYdim below.
%
% ---NOTE on a validated pathology and how it's handled: for series that
% sit *exactly* on a low-dimensional manifold with no measurement noise
% (e.g. a textbook periodic sine wave, or the logistic map, both tested
% directly), embedding in m=3 dimensions over-embeds by one or two
% spurious directions, and the local-Jacobian estimation becomes
% ill-conditioned -- producing large, spurious *positive* exponents that
% look like chaos but aren't (e.g. a clean sine wave gave LE1 = 0.46 at
% m=3, vs. the correct near-zero pair obtained at m=2). This is very
% unlikely to affect real (noisy) empirical data, but does bite on clean
% synthetic benchmark series (of which hctsa's own test suites include a
% few, e.g. dynamical-system-generated series). Since even a tiny amount
% of noise resolves it (0.1% of the series' own std was already enough
% to fix the sine-wave case above, giving LE1 = -0.0003), a fixed,
% reproducibly-seeded dither at that level is added internally before
% embedding -- negligible next to any real measurement noise, but enough
% to break the exact degeneracy.
%
% ---INPUTS:
% y, the input time series
%
% tauMethod, the time-delay for the embedding: an integer, or 'ac'/'mi'
%       for the first zero-crossing of the autocorrelation function or
%       first minimum of the automutual information (default: 'mi', as
%       used elsewhere for delay-embedding-based operations, e.g.
%       NL_EmbedPCA)
%
% m, the embedding dimension, and hence the number of Lyapunov exponents
%       estimated (default: 3 -- the minimum dimension in which a flow
%       can be chaotic at all (Poincare-Bendixson), and, in practice,
%       about the most that can be reliably resolved from typical
%       finite real-world data regardless of the true attractor's
%       dimension -- higher m both costs more and, per the note above,
%       becomes progressively less reliable to estimate). Fixed (rather
%       than auto-determined per series, e.g. via false-nearest-
%       neighbours, as NL_LargestLyap does for its own single-number
%       output) so that LE1/LE2/LE3 mean the same thing, and are
%       directly comparable, across every time series -- hctsa's feature
%       matrix needs a fixed number of output columns.
%
% kNN, the number of nearest neighbors used to estimate each local
%       Jacobian (default: 30, TISEAN's own default)
%
% maxN, the maximum number of samples to consider. Time series longer
%       than this are cropped to their first maxN points (default:
%       10000, cf. NL_RQA/NL_ZeroOneTest's maxN, same rationale -- cost
%       scales with series length, roughly linearly to slightly
%       superlinearly in practice; cost does not appreciably depend on
%       m). Set to 'full' to disable cropping.
%
% ---OUTPUTS:
% LE1, LE2, LE3, the three Lyapunov exponents, in descending order (LE1
%       is the largest/most expanding direction; LE3, typically the most
%       negative, is usually also the least reliably estimated of the
%       three)
% sumPos, the sum of the positive exponents -- an upper bound on the
%       Kolmogorov-Sinai entropy via Pesin's identity
% numPos, the number of positive exponents
% sumAll, the sum of LE1+LE2+LE3 -- the average phase-space
%       expansion/contraction rate (negative for a dissipative system;
%       for a flow this relates to the trace of the system's Jacobian)
% KYdim, the Kaplan-Yorke dimension: with exponents sorted descending
%       and k the largest number of leading exponents that still sum to
%       >= 0, KYdim = k + (sum of the first k exponents)/|LE(k+1)| -- a
%       conjectured estimate of the attractor's (fractal, box-counting)
%       dimension from the exponent spectrum alone. Capped at 3 if even
%       LE1+LE2+LE3 is >= 0 (formula doesn't apply; indicates m=3 wasn't
%       enough to see the attractor's contracting directions at all)

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
if nargin < 2 || isempty(tauMethod)
    tauMethod = 'mi';
end
if nargin < 3 || isempty(m)
    m = 3;
end
if nargin < 4 || isempty(kNN)
    kNN = 30;
end
if nargin < 5 || isempty(maxN)
    maxN = 10000;
end

y = y(:);
N = length(y);
if ischar(maxN) && strcmp(maxN, 'full')
    slowThreshold = 50000;
    if N > slowThreshold
        warning('Time series (%u samples) exceeds %u with maxN=''full''; computation may be slow', N, slowThreshold);
    end
elseif N > maxN
    warning('Time series (%u samples) exceeds maxN = %u; analyzing the first %u samples', N, maxN, maxN);
    y = y(1:maxN);
    N = maxN;
end

% ------------------------------------------------------------------------------
%% Tiny reproducible dither, to avoid the exact-low-dimensional-manifold
%% pathology described above (negligible relative to real measurement noise;
%% locally seeded and restored so this has no effect on the caller's rng state
%% and the dither itself is identical, hence reproducible, across calls)
% ------------------------------------------------------------------------------
rngState = rng(42, 'twister');
y = y + 0.001 * std(y) * randn(size(y));
rng(rngState);

% ------------------------------------------------------------------------------
%% Resolve the embedding delay tau (m is already fixed/numeric)
% ------------------------------------------------------------------------------
tm = BF_Embed(y, tauMethod, m, true);
tau = tm(1);
if isnan(tau)
    warning('Could not determine a suitable time delay for this time series');
    out = NaN; return
end

% ------------------------------------------------------------------------------
%% Run the TISEAN code, lyap_spec
% ------------------------------------------------------------------------------
filePath = BF_WriteTempFile(y);
outFilePath = [filePath '.lyaps'];

[~, res] = BF_TiseanSystem(sprintf('lyap_spec -m1,%u -d%u -k%u -o %s %s', ...
                          m, tau, kNN, outFilePath, filePath));

if exist(filePath, 'file'), delete(filePath); end

if isempty(res) || ~isempty(regexp(res, 'command not found', 'once'))
    if exist(outFilePath, 'file'), delete(outFilePath); end
    error('Call to TISEAN function ''lyap_spec'' failed.');
end

if ~exist(outFilePath, 'file')
    error('TISEAN function ''lyap_spec'' did not produce a .lyaps output file.');
end

fileInfo = dir(outFilePath);
if fileInfo.bytes == 0
    delete(outFilePath);
    out = NaN; return
end

% ------------------------------------------------------------------------------
%% Parse the output: a single data row (the converged exponent estimates),
%% plus a handful of '#'-prefixed comment lines with diagnostics we don't need
%% (TISEAN's own KYdim estimate is recomputed independently below instead of
%% parsed, for full control over edge cases like k==m)
% ------------------------------------------------------------------------------
fid = fopen(outFilePath, 'r');
dataLine = '';
while true
    tline = fgetl(fid);
    if ~ischar(tline), break; end
    if ~isempty(tline) && tline(1) ~= '#'
        dataLine = tline;
        break
    end
end
fclose(fid);
delete(outFilePath);

if isempty(dataLine)
    out = NaN; return
end

vals = sscanf(dataLine, '%f');
if length(vals) ~= m + 1 % n_used, then m exponents
    out = NaN; return
end
LE = vals(2:end); % already descending (largest/most expanding first)

if any(isnan(LE)) || any(isinf(LE))
    out = NaN; return
end

% ------------------------------------------------------------------------------
%% Output statistics
% ------------------------------------------------------------------------------
out.LE1 = LE(1);
out.LE2 = LE(2);
out.LE3 = LE(3);

out.sumPos = sum(LE(LE > 0));
out.numPos = sum(LE > 0);
out.sumAll = sum(LE);

cumLE = cumsum(LE);
posK = find(cumLE < 0, 1, 'first') - 1;
if isempty(posK) % even the full cumulative sum stays >= 0
    out.KYdim = m;
elseif posK == 0
    out.KYdim = 0;
else
    out.KYdim = posK + cumLE(posK) / abs(LE(posK + 1));
end

end
