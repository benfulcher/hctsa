function out = CO_Quantilogram(y, lag)
% CO_Quantilogram   Serial dependence of the series being below its own quantiles.
%
% Computes the quantilogram of Linton and Whang (2007): the autocorrelation, at a
% given lag, of the 'hit' process that marks when the series is below its
% alpha-quantile, for a range of quantile levels alpha from the lower tail
% through the center to the upper tail. At level alpha, the hit process is
% h_t = 1(y_t <= q_alpha), where q_alpha is the sample alpha-quantile, and the
% quantilogram is the sample autocorrelation of h_t at the lag:
% sum_t (h_t - hbar)(h_{t+lag} - hbar) / sum_t (h_t - hbar)^2.
% It is positive when excursions below the alpha-quantile cluster in time and
% negative when they alternate with excursions above it. For an independent
% series, it is zero at every level, with a standard error of about 1/sqrt(N).
% Looking at the lower and upper quantile levels separately captures dependence
% in the tails (volatility clustering, or bursts of extreme values that follow
% one another), at the center (directional persistence), and any asymmetry
% between them. Because the series enters only through whether each value lies
% below a quantile, the result at a fixed lag is unchanged by any increasing
% monotonic rescaling of the series, and is not affected by outliers (a lag set
% by the timescale of the series is found from the series itself, so it can
% change with the rescaling).
%
% ---INPUTS:
% y, the input time series
% lag, the time lag (default: 1). Can be a positive integer, 'ac' (the first
%       zero-crossing of the autocorrelation function of the series), 'ac1e'
%       (the floor of its first 1/e crossing), or 'mi' (the smaller of the first
%       minimum of the Kraskov automutual information and the 'ac1e' delay): lags
%       set by the timescale of the series (see BF_GetTau).
%
% ---OUTPUTS:
% A structure with the quantilogram at the lag for each quantile level, in
% fields q05, q10, q25, q50, q75, q90 and q95, for alpha = 0.05, 0.10, 0.25,
% 0.50, 0.75, 0.90 and 0.95, respectively. A level for which the hit process is
% constant (a constant series) or the lag is not smaller than N/2 is NaN.
%
% ---REFERENCES:
% O. Linton and Y.-J. Whang, "The quantilogram: With an application to evaluating
% directional predictability", J. Econometrics 141, 250 (2007).
%
% ---NOTES:
% The sample quantile is the round(alpha*N)-th smallest value, and the hit
% process is centered at its sample mean rather than at alpha, which is the same
% for a continuous-valued series (they differ by at most 1/N) and stays
% well-defined for series with tied values.
%
% The central quantile levels measure the same persistence as the quantile-state
% transition probabilities and the up-down motif frequencies (SB_TransitionMatrix,
% SB_MotifTwo): at the adaptive lag set by 'mi' the hit processes for q25, q50 and
% q75 are transforms of SB_TransitionMatrix stay-probabilities (maximum absolute
% Spearman correlation 0.96-0.97 on real series), so hctsa does not register them.
% It registers only the tail levels: q05, q10, q90 and q95 at the adaptive lag set
% by 'mi' (maximum absolute correlation with the SB_TransitionMatrix features 0.88),
% which grows with the series' own timescale and so is less dependent on the
% sampling rate (a series that decorrelates within a step, such as a chaotic map,
% keeps a lag of 1), and q05 and q95 at lag 1 (correlation 0.93 with the
% SB_TransitionMatrix_51 diagonal entries, accepted as interpretable). All seven
% levels remain computed and available.

% ------------------------------------------------------------------------------
% Copyright (C) 2026, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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
if nargin < 2 || isempty(lag)
    lag = 1;
end

alphas = [0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95];
fieldNames = {'q05','q10','q25','q50','q75','q90','q95'};
out = cell2struct(num2cell(nan(1, numel(alphas))), fieldNames, 2);

y = y(:);
y = y(isfinite(y));
N = length(y);

if ischar(lag) || isstring(lag)
    switch lag
    case 'ac'
        lag = CO_FirstCrossing(y, 'ac', 0, 'discrete');
    case {'ac1e', 'mi'}
        lag = BF_GetTau(y, lag); % adaptive delay: see BF_GetTau
    otherwise
        error('Unknown lag option ''%s'': use a positive integer, ''ac'', ''ac1e'', or ''mi''', lag)
    end
end
if isnan(lag) || lag < 1 || lag >= N/2
    return % no lag defined, or too few pairs of observations
end

% ------------------------------------------------------------------------------
%% Autocorrelation of the quantile-hit process at each quantile level
% ------------------------------------------------------------------------------
ySorted = sort(y);
for i = 1:numel(alphas)
    q = ySorted(max(1, round(alphas(i)*N))); % sample alpha-quantile
    h = double(y <= q);
    h = h - mean(h);
    denom = sum(h.^2);
    if denom > 0
        out.(fieldNames{i}) = sum(h(1:N-lag).*h(1+lag:N))/denom;
    end
end

end
