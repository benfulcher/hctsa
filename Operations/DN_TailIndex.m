function out = DN_TailIndex(y, tailFrac)
% DN_TailIndex   Tail index of the distribution of values: how heavy its tails are.
%
% Estimates the extreme-value index of the marginal distribution of the values,
% ignoring their temporal order, from the most extreme tailFrac proportion of
% observations at each end. The extreme-value index xi describes how the
% probability of very large values decays: xi > 0 is a heavy (power-law) tail
% with tail index 1/xi (the Student-t distribution with nu degrees of freedom
% has xi = 1/nu), xi = 0 an exponential-type tail (including the Gaussian), and
% xi < 0 a tail with a finite upper end point (a uniform distribution has
% xi = -1). Unlike kurtosis, the estimates depend only on the largest few
% percent of values and are not driven by a single outlier.
%
% Tails are defined relative to the median, so that the Hill estimator, which
% needs positive values, applies to any series: the upper tail is the values
% above the median (as distances above it), the lower tail is the values below
% the median (as distances below it), and the absolute tail (used by the moment
% estimator) is the distance of every value from the median, |y - median(y)|.
% Each estimator uses the k largest tail values, where k = round(tailFrac*N) is
% a fixed proportion of the series length N (so estimates are not trivially
% length-dependent), and takes logarithms of ratios to the (k+1)-th largest,
% which makes it unchanged by any rescaling of the series.
%
% Estimators:
% (i) the Hill (1975) estimator, the mean log-ratio of the k largest tail values
%       to the (k+1)-th largest. It assumes a heavy tail (xi > 0) and is biased
%       upward in samples from a light-tailed distribution (about 0.2 for a Gaussian
%       at tailFrac = 0.05, even though xi = 0, and always positive), and also
%       for skewed distributions because tails are measured from the median;
% (ii) the moment estimator of Dekkers, Einmahl and de Haan (1989), which
%       extends the Hill estimator to any sign of xi, including light-tailed
%       distributions with a finite end point (xi < 0);
% (iii) the shape parameter of a generalized Pareto distribution fitted by
%       maximum likelihood to the k exceedances over the (k+1)-th largest value
%       (the peaks-over-threshold method; Pickands, 1975), using gpfit from
%       MATLAB's Statistics Toolbox.
%
% ---INPUTS:
% y, the input time series
% tailFrac, the proportion of the N observations in the tail used for estimation
%       (default: 0.05, i.e., the largest 5% of values; the generalized Pareto
%       threshold is then the 95th percentile of the series). Needs at least 10
%       tail values (k >= 10), so series shorter than 10/tailFrac samples give NaN.
%
% ---OUTPUTS:
% A structure with fields:
% hillUpper, Hill estimate of xi for the upper tail (values above the median)
% hillLower, Hill estimate of xi for the lower tail (values below the median)
% hillAsym, tail asymmetry: hillUpper - hillLower (positive when the upper tail
%       is heavier than the lower tail)
% momentAbs, Dekkers-Einmahl-de Haan moment estimate of xi for the absolute tail
% gpdUpper, generalized Pareto shape parameter fitted to the upper tail
% gpdLower, generalized Pareto shape parameter fitted to the lower tail
% gpdAsym, tail asymmetry: gpdUpper - gpdLower
% Any quantity that cannot be computed (too few tail values, a tie between the
% k-th and (k+1)-th largest tail values, as happens for strongly discretized
% series, or a failed fit) is returned as NaN.
%
% ---REFERENCES:
% B.M. Hill, "A simple general approach to inference about the tail of a
% distribution", Ann. Stat. 3, 1163 (1975).
% A.L.M. Dekkers, J.H.J. Einmahl, L. de Haan, "A moment estimator for the index
% of an extreme-value distribution", Ann. Stat. 17, 1833 (1989).
% J. Pickands III, "Statistical inference using extreme order statistics",
% Ann. Stat. 3, 119 (1975).
%
% ---NOTES:
% With k tail values, the standard error of the Hill estimate is about xi/sqrt(k),
% and that of the moment and generalized Pareto estimates is about 0.2 for
% k = 50 (N = 1000 at tailFrac = 0.05), whatever xi is. All three are biased by the
% choice of k when the distribution is not exactly Pareto-tailed. The values are
% best read as a relative ranking of tail heaviness between series rather than as
% precise estimates of xi. In samples of N = 1000 from a Student-t distribution
% with 3 to 10 degrees of freedom (xi = 1/nu), the Hill estimate reads high by
% 0.05 to 0.15, and the moment and generalized Pareto estimates read low by about
% 0.1; for nu <= 2 the three agree to within about 0.1.
%
% hctsa registers hillUpper for tailFrac = 0.05 but not for 0.10: that field is nearly
% redundant with an outlier statistic (Spearman correlation 0.96 across series with
% DN_RemovePoints_max_01_saturate_mean).

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
% (gpfit, used for the generalized Pareto fits, is in the Statistics and Machine
% Learning Toolbox)
BF_CheckToolbox('statistics_toolbox');

if nargin < 2 || isempty(tailFrac)
    tailFrac = 0.05;
end

fieldNames = {'hillUpper','hillLower','hillAsym','momentAbs', ...
                'gpdUpper','gpdLower','gpdAsym'};
out = cell2struct(num2cell(nan(1, numel(fieldNames))), fieldNames, 2);

y = y(:);
y = y(isfinite(y));
N = length(y);
k = round(tailFrac*N); % number of values in each tail
if k < 10 || k >= floor(N/2)
    return % too few tail values to estimate a tail index
end

% ------------------------------------------------------------------------------
%% Tails relative to the median (in descending order)
% ------------------------------------------------------------------------------
dev = y - median(y);
tailUp = sort(dev(dev > 0), 'descend'); % distances above the median
tailLo = sort(-dev(dev < 0), 'descend'); % distances below the median
tailAbs = sort(abs(dev), 'descend'); % distances from the median

% ------------------------------------------------------------------------------
%% Hill and moment estimators
% ------------------------------------------------------------------------------
out.hillUpper = hillEstimate(tailUp, k);
out.hillLower = hillEstimate(tailLo, k);
out.hillAsym = out.hillUpper - out.hillLower;
out.momentAbs = momentEstimate(tailAbs, k);

% ------------------------------------------------------------------------------
%% Generalized Pareto fits to exceedances over the (k+1)-th largest value
% ------------------------------------------------------------------------------
ys = sort(y, 'descend');
out.gpdUpper = gpdShape(ys(1:k) - ys(k + 1)); % upper exceedances
ys = sort(y, 'ascend');
out.gpdLower = gpdShape(ys(k + 1) - ys(1:k)); % lower exceedances
out.gpdAsym = out.gpdUpper - out.gpdLower;

end

% ------------------------------------------------------------------------------
function xi = hillEstimate(s, k)
    % s: positive tail values in descending order
    if numel(s) <= k || s(k + 1) <= 0 || s(k) == s(k + 1)
        xi = NaN; % too few values, or a tie at the threshold
    else
        xi = mean(log(s(1:k))) - log(s(k + 1));
    end
end

% ------------------------------------------------------------------------------
function xi = momentEstimate(s, k)
    % Dekkers-Einmahl-de Haan moment estimator
    if numel(s) <= k || s(k + 1) <= 0 || s(k) == s(k + 1)
        xi = NaN; return % too few values, or a tie at the threshold
    end
    logExcess = log(s(1:k)) - log(s(k + 1));
    M1 = mean(logExcess);
    M2 = mean(logExcess.^2);
    if M2 <= 0
        xi = NaN; return
    end
    xi = M1 + 1 - 0.5/(1 - M1^2/M2);
    if ~isfinite(xi)
        xi = NaN;
    end
end

% ------------------------------------------------------------------------------
function xi = gpdShape(exceed)
    % Maximum-likelihood shape of a generalized Pareto fit (Statistics Toolbox)
    if any(exceed <= 0)
        xi = NaN; return % ties between the tail values and the threshold
    end
    warnState = warning('off', 'all');
    try
        parmhat = gpfit(exceed);
        xi = parmhat(1);
    catch
        xi = NaN;
    end
    warning(warnState);
end
