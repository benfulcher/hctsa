function out = CP_WaveletVarChg(y, wName, level, maxnchpts, minDelay)
% CP_WaveletVarChg   Number of variance change points found from a wavelet detail signal.
%
% The series is decomposed with a mother wavelet to a given level (Wavelet
% Toolbox), and the detail signal at that level is reconstructed. The 2% of its
% values with the largest magnitude are replaced by the mean of the detail
% signal, so that single large events do not dominate. For k = 0, ..., maxnchpts,
% the best segmentation of the (centered) detail signal into k + 1 segments of
% constant variance, each at least minDelay samples long, is found exactly by
% dynamic programming on the Gaussian variance contrast of Lavielle (1999),
% J_k = sum_j n_j*log(S_j/n_j), where segment j has n_j samples and sum of
% squares S_j. The number of change points is then the k that minimizes the
% modified BIC of Zhang and Siegmund (2007),
%       J_k/tau + sum_j log(n_j) + k*log(N),
% in which the contrast is divided by tau, the variance inflation of segment
% mean squares caused by the detail signal's own serial dependence and heavy
% tails (see NOTES). The output is that number.
% Needs the Wavelet Toolbox. Returns NaN if level is too large for the series.
%
% ---INPUTS:
% y, the input time series
%
% wName, the name of the mother wavelet to analyze the data with: e.g., 'db3',
%       'sym2', cf. Wavelet Toolbox Documentation for details (default: 'db3')
%
% level, the level of wavelet decomposition (default: 3); 'max' uses the maximum
%       level for the series length
%
% maxnchpts, the maximum number of change points (default: 5)
%
% minDelay, the minimum delay between consecutive change points, in samples; a
%       value between 0 and 1 is taken as a proportion of the time-series
%       length, e.g., 0.02 ensures that change points are separated by at least
%       2% of the time-series length (default: 0.01)
%
% ---OUTPUTS:
% a scalar: the number of change points, from 0 to maxnchpts, or NaN if level is
%       too large for wavelet wName on this series
%
% ---NOTES:
% The detail signal is band-limited, so it is strongly autocorrelated even when
% the input is white noise, and an information criterion that treats its samples
% as independent reports change points in nearly every stationary series. For a
% stationary Gaussian signal with autocorrelation rho(h), the mean square of a
% segment of n samples has a variance inflated by the factor
% 1 + 2*sum_h rho(h)^2 relative to n independent samples. tau estimates this
% factor from the detail signal itself, with a Bartlett taper up to lag
% L = min(200, floor(N/4)),
%       tau = 1 + 2*sum_{h=1..L} (1 - h/(L+1))*rho(h)^2 (at least 1),
% multiplied by max(1, (kurtosis - 1)/2) for heavy tails (kurtosis 3 for
% Gaussian data). Dividing the contrast by tau counts N/tau effective
% observations, so that stationary series (including autocorrelated and
% heavy-tailed ones) only rarely yield a spurious change point.
% The deeper the decomposition level, the fewer effective observations the
% detail signal has (tau is about 2, 4, 8 and 15 at levels 2-5 for white noise),
% so only larger variance changes are detected there.
%
% ---REFERENCES:
% M. Lavielle, "Detection of multiple changes in a sequence of dependent
% variables", Stochastic Processes and their Applications 83(1), 79 (1999).
% DOI: 10.1016/S0304-4149(99)00023-X
% N. R. Zhang and D. O. Siegmund, "A modified Bayes information criterion with
% applications to the analysis of comparative genomic hybridization data",
% Biometrics 63(1), 22 (2007). DOI: 10.1111/j.1541-0420.2006.00662.x

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

% Check that a Wavelet Toolbox license is available:
BF_CheckToolbox('wavelet_toolbox');

% ------------------------------------------------------------------------------
%% Check Inputs
% ------------------------------------------------------------------------------
N = length(y); % time-series length

if nargin < 2 || isempty(wName)
	wName = 'db3'; % default wavelet
end

if nargin < 3 || isempty(level)
	level = 3; % level of wavelet decomposition
end
if strcmp(level, 'max')
	level = wmaxlev(N, wName);
end

if nargin < 4 || isempty(maxnchpts)
	maxnchpts = 5; % maximum number of change points
end

if nargin < 5 || isempty(minDelay)
	minDelay = 0.01; % 1% of the time series length
end
if (minDelay > 0) && (minDelay < 1)
	minDelay = ceil(minDelay * N);
end

if wmaxlev(N, wName) < level
	warning('Chosen level, %u, is too large for wavelet %s on a series of length %u', level, wName, N);
	out = NaN; return
end

% ------------------------------------------------------------------------------
% 1. Detail signal at the chosen level
% ------------------------------------------------------------------------------
[c, l] = wavedec(y, level, wName);
det = wrcoef('d', c, l, wName, level);

% ------------------------------------------------------------------------------
% 2. Replace 2% of the greatest (absolute) values by the mean
% ------------------------------------------------------------------------------
x = sort(abs(det));
v2p100 = x(fix(length(x) * 0.98));
det(abs(det) > v2p100) = mean(det);

% ------------------------------------------------------------------------------
% 3. Best segmentations with 0..maxnchpts change points, and the number to keep
% ------------------------------------------------------------------------------
det = det(:);
if numel(det) < 2 || ~all(isfinite(det)) || maxnchpts < 1
	out = 0; return
end
[J, segEnds] = VarianceSegmentations(det, maxnchpts, minDelay);
out = NumChangesMBIC(J, segEnds, numel(det), VarianceInflation(det));

end
% ------------------------------------------------------------------------------
function [J, segEnds] = VarianceSegmentations(x, K, d)
	% Exact minimum of the Gaussian variance contrast, sum_j n_j*log(S_j/n_j), over
	% segmentations of x into k + 1 segments of at least d samples (k = 0..K), by
	% dynamic programming. J(k+1) is the minimum (Inf if infeasible); segEnds(k+1,1:k+1)
	% are the segment ends (the last is N; zeros if infeasible).
	N = numel(x);
	d = max(floor(d), 1);
	x = x - mean(x);
	cs = [0; cumsum(x.^2)];
	% F(m+1, e+1): minimal contrast of the first e samples in m segments
	F = Inf(K + 2, N + 1);
	F(1, 1) = 0;
	back = zeros(K + 2, N + 1);
	for e = d:N
		sMax = e - d; % (the last segment has at least d samples)
		s = 0:sMax;
		n = e - s;
		c = n .* log(max(cs(e + 1) - cs(s + 1)', 0) ./ n);
		cand = F(1:K + 1, 1:sMax + 1) + c;
		if any(c == -Inf) % (a zero-variance segment: Inf + -Inf counts as infeasible)
			cand(isnan(cand)) = Inf;
		end
		[F(2:end, e + 1), idx] = min(cand, [], 2);
		back(2:end, e + 1) = idx - 1;
	end
	J = Inf(K + 1, 1);
	segEnds = zeros(K + 1, K + 1);
	for k = 0:K
		J(k + 1) = F(k + 2, N + 1);
		if ~isfinite(J(k + 1))
			continue
		end
		e = N;
		ends = zeros(1, k + 1);
		for mm = k + 1:-1:1
			ends(mm) = e;
			e = back(mm + 1, e + 1);
		end
		segEnds(k + 1, 1:k + 1) = ends;
	end
end
% ------------------------------------------------------------------------------
function tau = VarianceInflation(x)
	% Variance inflation of segment mean squares from serial dependence (Bartlett-tapered
	% sum of squared autocorrelations) and heavy tails; at least 1.
	N = numel(x);
	xc = x - mean(x);
	v = mean(xc.^2);
	if N < 8 || ~(v > 0)
		tau = 1; return
	end
	L = min(200, floor(N/4));
	f = fft(xc, 2*N);
	a = real(ifft(f .* conj(f)));
	rho = a(2:L + 1) / a(1);
	h = (1:L)';
	tau = 1 + 2*sum((1 - h/(L + 1)) .* rho.^2);
	kurt = mean(xc.^4) / v^2;
	tau = max(tau, 1) * max(1, (kurt - 1)/2);
end
% ------------------------------------------------------------------------------
function kBest = NumChangesMBIC(J, segEnds, N, tau)
	% k minimizing J_k/tau + sum_j log(n_j) + k*log(N) (Zhang and Siegmund, 2007)
	best = Inf;
	kBest = 0;
	for k = 0:numel(J) - 1
		if ~isfinite(J(k + 1))
			break
		end
		n = diff([0, segEnds(k + 1, 1:k + 1)]);
		if any(n <= 0)
			continue
		end
		obj = J(k + 1)/tau + sum(log(n)) + k*log(N);
		if obj < best
			best = obj;
			kBest = k;
		end
	end
end
