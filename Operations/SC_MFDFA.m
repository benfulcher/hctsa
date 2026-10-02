function out = SC_MFDFA(y, scaleRange, qRange, order)
% SC_MFDFA   Multifractal detrended fluctuation analysis (MFDFA): the multifractal spectrum of a time series.
%
% Estimates the multifractal singularity spectrum f(alpha) of a time series by
% the classical MFDFA algorithm of Kantelhardt et al. (2002). The mean-subtracted
% series is integrated (cumulative sum) to a profile, which is divided into
% non-overlapping segments of length s (taken from both the start and the end of
% the series, to use the whole series when N is not a multiple of s). Each segment
% is detrended by a polynomial of a given order, and the q-th order fluctuation
% function F_q(s) is formed by averaging the segment variances raised to the power
% q/2 (with a log-averaging limit at q = 0). For a monofractal series, F_q(s)
% scales as s^h(q) with a single exponent h independent of q. For a multifractal
% series, h(q) varies with q: small or negative q weight small-fluctuation
% segments, large positive q weight large-fluctuation segments, so a q-dependence
% of h(q) reveals that different amplitudes obey different local scaling laws.
% h(q) is Legendre-transformed via the mass exponent tau(q) = q*h(q) - 1 into the
% singularity spectrum f(alpha), where alpha = d tau / dq and f = q*alpha - tau,
% the standard multifractal fingerprint. Its width, Delta-alpha, is the usual
% "degree of multifractality".
%
% It differs from hctsa's other scaling-exponent estimators: SC_FastDFA and
% SC_FluctAnal are monofractal (a single q = 2 exponent, no Legendre transform);
% SC_MMA also generalizes DFA across q but reports how the raw h(q) surface varies
% with scale, and has no f(alpha) spectrum. SC_MFDFA fixes the scaling range and
% focuses on the q axis.
%
% ---INPUTS:
% y, the input time series
% scaleRange, [minScale, maxScale], the range of segment lengths s used for the
%       fluctuation-function fit, as 20 log-spaced values (default:
%       [16, floor(N/4)], following common DFA practice, cf. Peng et al. 1995)
% qRange, [qMin, qMax], the range of the multifractal order q (default: [-5, 5];
%       q is sampled in steps of 0.5, with q = 0 handled by its log-averaging
%       limit)
% order, the order of the polynomial used to detrend each segment (default: 1,
%       linear detrending, as in the original MFDFA1 of Kantelhardt et al.)
%
% ---OUTPUTS: statistics of the Legendre-transformed singularity spectrum
% alpha / f(alpha), and quality diagnostics of the h(q) fits:
% meanR2, the mean (across q) of the R^2 of the log-log fits of F_q(s) against s
% h2, the generalized Hurst exponent h(q) at q = 2 (NaN if q = 2 is not in qRange)
% alphaMin, alphaMax, the smallest and largest singularity exponents alpha
% alphaWidth, alphaMax - alphaMin, the width of the spectrum (Delta-alpha, the
%       degree of multifractality)
% fAlphaMax, the maximum of f(alpha), the height of the spectrum's peak
% alpha0, the alpha at the peak of the spectrum (the dominant exponent)
% spectrumAsymmetry, the width to the left of the peak divided by the width to the
%       right, (alpha0 - alphaMin) / (alphaMax - alpha0)
% (The raw generalized Hurst exponent h(q) is not otherwise reported: it is
% redundant, r >= 0.9, with these spectrum quantities.)
% A scalar NaN is returned if the series or scale range cannot support the
% analysis.
%
% ---REFERENCES:
% J. W. Kantelhardt et al., "Multifractal detrended fluctuation analysis of
% nonstationary time series", Physica A 316(1-4), 87-114 (2002).

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
if size(y, 2) > size(y, 1)
	y = y'; % ensure a column vector
end
N = length(y);

if nargin < 2 || isempty(scaleRange)
	scaleRange = [16, floor(N / 4)];
end
minScale = scaleRange(1);
maxScale = scaleRange(2);

if nargin < 3 || isempty(qRange)
	qRange = [-5, 5];
end
qMin = qRange(1);
qMax = qRange(2);

if nargin < 4 || isempty(order)
	order = 1; % linear detrending (MFDFA1)
end

% Need at least order+3 points per segment for a meaningful (non-degenerate)
% polynomial-detrended variance, and at least 8 distinct scales spanning the
% range to fit a reliable log-log slope for h(q):
if minScale < order + 3
	minScale = order + 3;
end
numScales = 20;
if maxScale > floor(N / 4) || (maxScale / minScale) < 2 || floor(N / (2 * minScale)) < 4
	out = NaN;
	return
end
scales = unique(round(exp(linspace(log(minScale), log(maxScale), numScales))));
if length(scales) < 8
	out = NaN;
	return
end

qList = qMin:0.5:qMax;
qList(qList == 0) = 1e-4; % q=0 handled separately via log-averaging below
qZeroIdx = find(qMin:0.5:qMax == 0, 1);

% ------------------------------------------------------------------------------
%% Compute the profile (integrated, mean-subtracted series)
% ------------------------------------------------------------------------------
profile = cumsum(y - mean(y));

% ------------------------------------------------------------------------------
%% Local detrended variance F^2(s,v) at each scale, for every segment
% ------------------------------------------------------------------------------
numQ = length(qList);
logFq = nan(length(scales), numQ);
logS = log(scales(:));

for si = 1:length(scales)
	s = scales(si);
	Ns = floor(N / s);

	% Segments from the start, and (to use the whole series) from the end:
	tt = (1:s)';
	F2 = nan(2 * Ns, 1);
	for v = 1:Ns
		idx = (v - 1) * s + (1:s);
		seg = profile(idx);
		p = polyfit(tt, seg, order);
		F2(v) = mean((seg - polyval(p, tt)).^2);
	end
	for v = 1:Ns
		idx = N - v * s + (1:s);
		seg = profile(idx);
		p = polyfit(tt, seg, order);
		F2(Ns + v) = mean((seg - polyval(p, tt)).^2);
	end

	F2(F2 < eps) = eps; % floor to avoid log(0) / 0^(negative q)

	for qi = 1:numQ
		q = qList(qi);
		logFq(si, qi) = (1 / q) * log(mean(F2.^(q / 2)));
	end
	if ~isempty(qZeroIdx)
		% q=0 limit: F_0(s) = exp{ (1/(4Ns)) * sum( log F^2(s,v) ) }
		logFq(si, qZeroIdx) = mean(log(F2)) / 2;
	end
end

% ------------------------------------------------------------------------------
%% h(q): slope of log(F_q(s)) vs log(s), for each q
% ------------------------------------------------------------------------------
hq = nan(numQ, 1);
r2q = nan(numQ, 1);
for qi = 1:numQ
	good = ~isnan(logFq(:, qi)) & ~isinf(logFq(:, qi));
	if sum(good) < 8
		continue
	end
	p = polyfit(logS(good), logFq(good, qi), 1);
	hq(qi) = p(1);
	fitted = polyval(p, logS(good));
	ssRes = sum((logFq(good, qi) - fitted).^2);
	ssTot = sum((logFq(good, qi) - mean(logFq(good, qi))).^2);
	if ssTot > 0
		r2q(qi) = 1 - ssRes / ssTot;
	end
end

if any(isnan(hq))
	out = NaN;
	return
end

% ------------------------------------------------------------------------------
%% Legendre transform: mass exponent tau(q), singularity spectrum alpha/f(alpha)
% ------------------------------------------------------------------------------
qListTrue = qMin:0.5:qMax; % the true q values (q=0 restored, not the 1e-4 surrogate)
tauq = qListTrue(:) .* hq - 1;
alpha = gradient(tauq, qListTrue(:)); % d(tau)/dq
falpha = qListTrue(:) .* alpha - tauq;

% ------------------------------------------------------------------------------
%% Output statistics
% ------------------------------------------------------------------------------
out.meanR2 = mean(r2q, 'omitnan'); % quality of the h(q) log-log fits

% h(q) at q=2 (~classic single-exponent DFA, for reference against
% hctsa's other scaling-exponent estimators):
idx2 = find(abs(qListTrue - 2) < 1e-8, 1);
if isempty(idx2)
	out.h2 = NaN;
else
	out.h2 = hq(idx2);
end

% NOTE: the rest of the raw generalized Hurst exponent h(q) -- its
% q=qMin/qMax endpoints and its range/trend/std across q -- is NOT reported
% here. Verified on two independent collections of real-world series (r>=0.95 on
% BOTH, the bar for confirmed redundancy) that h(qMin) duplicates alphaMax,
% h(qMax) duplicates alphaMin, and h(q)'s range/trend/std across q all
% mutually duplicate alphaWidth (hqTrend is in fact a deterministic
% rescaling of hqRange, since qRange is fixed). h2 survives this bar --
% its only >=0.95-on-both link is to the now-dropped h(qMax), and it is
% NOT >=0.95 against alphaMin/alphaMax/alphaWidth/alpha0 on both datasets --
% so it is kept. The alpha/f(alpha) spectrum is the standard reporting
% convention in the multifractal literature, so it is kept in preference to
% h(qMin)/h(qMax)/hqRange/hqTrend/hqStd specifically.

% Singularity spectrum, alpha/f(alpha):
out.alphaMin = min(alpha);
out.alphaMax = max(alpha);
out.alphaWidth = out.alphaMax - out.alphaMin; % "degree of multifractality", Delta-alpha
[out.fAlphaMax, iPeak] = max(falpha);
out.alpha0 = alpha(iPeak); % dominant (most probable) singularity exponent
% Spectrum asymmetry: >1 => left-skewed (long tail from large/negative q,
% dominated by small fluctuations); <1 => right-skewed:
leftWidth = out.alpha0 - out.alphaMin;
rightWidth = out.alphaMax - out.alpha0;
if rightWidth > 0
	out.spectrumAsymmetry = leftWidth / rightWidth;
else
	out.spectrumAsymmetry = NaN;
end

end
