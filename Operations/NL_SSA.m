function out = NL_SSA(y, L)
% NL_SSA  Singular Spectrum Analysis of a time series.
%
% Constructs the trajectory (Hankel) matrix of the time series using a
% window length L (i.e., a time-delay embedding with delay tau = 1), and
% performs an uncentered singular value decomposition of the result.
%
% Unlike NL_EmbedPCA (which centers the embedded data before decomposing it,
% and allows a general embedding delay), this implements classic "Basic SSA":
% a fixed delay of 1, no centering (so that a genuine trend is not removed
% before decomposition), and diagonal averaging ("Hankelization") of the
% leading elementary matrices back into component time series. Statistics
% are computed on the singular-value pairing structure, and on the
% reconstructed leading trend/oscillatory components themselves, rather than
% on the raw eigenvalue spectrum (which NL_EmbedPCA already covers).
%
% "Extracting qualitative dynamics from experimental data"
% D. S. Broomhead and G. P. King, Physica D 20(2-3) 217 (1986)
%
% "Analysis of Time Series Structure: SSA and Related Techniques"
% N. Golyandina, V. Nekrutkin, A. Zhigljavsky, Chapman & Hall/CRC (2001)
%
% ---INPUTS:
% y, the input time series
%
% L, the window length (default: floor(N/4)). Must satisfy 4 <= L <= floor(N/2).
%
% ---OUTPUTS: Statistics on the singular-value pairing/decay structure
% (gap1-gap5, sepIdx), and on the reconstructed leading block of components as
% a whole (trend strength trend_r2/trend_rho, dominant period of its most
% tightly-paired internal mode pairperiod, and w-correlation-based separability
% from the residual wcorr_leadresid). Individual eigentriples within a
% near-degenerate pair are not uniquely determined by the SVD, so all outputs
% are computed from basis-independent quantities: singular-value gaps, or sums
% of whole blocks of components rather than single components.

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
% Check inputs:
% ------------------------------------------------------------------------------
N = length(y);

if nargin < 2 || isempty(L)
    % Classic SSA guidance: L <= N/2; N/4 is a common default. But the SVD
    % below costs O(N*L^2), so letting L grow proportionally to N would make
    % this cubic in N for long time series -- cap the *default* window length
    % (in the same spirit as the maxL cap in NW_VisibilityGraph.m and the
    % maxN cap in NL_RQA.m). Explicitly-specified L is left uncapped: it is
    % the caller's choice (e.g., via an INP_mops_hctsa.txt variant).
    maxDefaultL = 200;
    L = min(floor(N/4), maxDefaultL);
end

if L < 4 || L > floor(N/2)
    fprintf(1, 'Window length L = %u is not suitable for a series of length %u\n', L, N);
    out = NaN; return
end

% ------------------------------------------------------------------------------
% Build the trajectory matrix (an embedding with delay tau = 1):
% ------------------------------------------------------------------------------
X = BF_Embed(y, 1, L, false); % K x L, K = N - L + 1

if isnan(X)
    fprintf(1, 'Could not construct a trajectory matrix for this time series\n');
    out = NaN; return
end

% ------------------------------------------------------------------------------
% Uncentered SVD of the trajectory matrix (do NOT mean-subtract: a genuine
% trend should survive into the leading component(s), unlike in NL_EmbedPCA):
% ------------------------------------------------------------------------------
[U, S, V] = svd(X, 'econ');
sigma = diag(S);
d = length(sigma);

if d < 4
    fprintf(1, 'Not enough singular values (%u) obtained for this window length\n', d);
    out = NaN; return
end

% ------------------------------------------------------------------------------
%% Singular-value pairing structure
% ------------------------------------------------------------------------------
% Small relative gaps between consecutive singular values are the classic SSA
% signature of a paired (oscillatory) mode; a large gap signals a
% well-separated, non-oscillatory mode (e.g., a trend), or the boundary
% between structured "signal" components and the noise continuum.
%
% NOTE: individual eigentriples within a near-degenerate pair (sigma_i approx
% sigma_i+1, the case that matters most, e.g., for a genuine oscillation) are
% NOT uniquely determined by the SVD -- any rotation within that subspace is
% an equally valid solution, so a single reconstructed component or its
% eigenvector is numerically arbitrary. Everything below therefore only ever
% uses gaps (basis-independent, from singular values alone) or SUMS of whole
% blocks of components (also basis-independent, since summing across a
% degenerate subspace is invariant to the rotation SVD happens to return).
maxCheck = min(6, d-1);
allGaps = (sigma(1:maxCheck) - sigma(2:maxCheck+1)) ./ sigma(1:maxCheck);
for i = 1:5
    if i <= maxCheck
        out.(sprintf('gap%u', i)) = allGaps(i);
    else
        out.(sprintf('gap%u', i)) = NaN;
    end
end

% ------------------------------------------------------------------------------
%% Leading structured block vs. residual
% ------------------------------------------------------------------------------
% sepIdx locates the largest relative drop in the spectrum: components
% 1:sepIdx are taken as the "structured" leading block, the rest as residual.
[~, sepIdx] = max(allGaps);
out.sepIdx = sepIdx;

leadElementary = zeros(size(X));
for i = 1:sepIdx
    leadElementary = leadElementary + sigma(i)*U(:,i)*V(:,i)';
end
c_lead = SSA_diagonalAverage(leadElementary, N);
resid = y(:) - c_lead; % exact, since diagonal-averaging the full X recovers y

% ------------------------------------------------------------------------------
%% Trend diagnostics on the leading block
% ------------------------------------------------------------------------------
t = (1:N)';
out.trend_rho = corr(c_lead, t, 'type', 'Spearman');

linFit = [t, ones(N,1)] \ c_lead;
c_leadhat = [t, ones(N,1)] * linFit;
ssRes = sum((c_lead - c_leadhat).^2);
ssTot = sum((c_lead - mean(c_lead)).^2);
if ssTot > 0
    out.trend_r2 = 1 - ssRes/ssTot;
else
    out.trend_r2 = NaN;
end

% ------------------------------------------------------------------------------
%% Period of the most tightly-paired mode within the leading block
% ------------------------------------------------------------------------------
if sepIdx >= 2
    [~, iStar] = min(allGaps(1:sepIdx-1));
    pairElementary = sigma(iStar)*U(:,iStar)*V(:,iStar)' + sigma(iStar+1)*U(:,iStar+1)*V(:,iStar+1)';
    c_pair = SSA_diagonalAverage(pairElementary, N);
    c_pair = c_pair - mean(c_pair);
    numSignChanges = sum(diff(sign(c_pair)) ~= 0);
    if numSignChanges > 0
        out.pairperiod = 2*(N-1)/numSignChanges;
    else
        out.pairperiod = NaN;
    end
else
    out.pairperiod = NaN;
end

% ------------------------------------------------------------------------------
%% w-correlation between the leading block and the residual (separability)
% ------------------------------------------------------------------------------
% |w-correlation| close to 0 means the leading block is cleanly separated from
% the rest of the series; close to 1 means the decomposition has not resolved
% a clean structured/residual split (e.g., for a series with no clear
% low-dimensional structure, such as white noise).
w = min([t, L*ones(N,1), N-t+1], [], 2);
wcorr = @(a,b) sum(w.*a.*b) / sqrt(sum(w.*a.^2) * sum(w.*b.^2));
out.wcorr_leadresid = wcorr(c_lead, resid);

% ------------------------------------------------------------------------------
    function c = SSA_diagonalAverage(Xi, N)
        % Reconstructs a length-N component series from an elementary matrix
        % by averaging over its anti-diagonals ("Hankelization").
        % (Local names deliberately avoid K/L: this is a nested function, so
        % reusing an outer-scope variable name here would silently overwrite it.)
        [numRows, numCols] = size(Xi);
        c = zeros(N, 1);
        for tt = 1:N
            jlo = max(1, tt-numRows+1);
            jhi = min(numCols, tt);
            jRange = jlo:jhi;
            iRange = tt - jRange + 1;
            c(tt) = mean(Xi(sub2ind([numRows, numCols], iRange, jRange)));
        end
    end

end
