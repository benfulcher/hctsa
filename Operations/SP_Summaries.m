function out = SP_Summaries(y, psdMeth, windowType, nf, doLogAbs)
% SP_Summaries  Statistics of the power spectrum of a time series
%
% The estimation can be done using a periodogram, using the periodogram code in
% Matlab's Signal Processing Toolbox, or a fast fourier transform, implemented
% using Matlab's fft code.
%
% ---INPUTS:
% y, the input time series
%
% psdMeth, the method of obtaining the spectrum from the signal:
%               (i) 'periodogram': periodogram
%               (ii) 'fft': fast fourier transform
%               (iii) 'welch': Welch's method
%
% windowType, the window to use:
%               (i) 'boxcar'
%               (ii) 'rect'
%               (iii) 'bartlett'
%               (iv) 'hann'
%               (v) 'hamming'
%               (vi) 'none'
%
% nf, the number of frequency components to include. If
%           empty (default), it's approximately length(y).
%
% doLogAbs, if 1, takes log amplitude of the signal before
%           transforming to the frequency domain.
%
% doPower, analyzes the power spectrum rather than amplitudes of a Fourier
%          transform
%
% ---OUTPUTS:
% Statistics summarizing various properties of the spectrum,
% including its maximum, minimum, spread, correlation, centroid, area in certain
% (normalized) frequency bands, moments of the spectrum, power-weighted
% moments of the frequency distribution, Shannon spectral
% entropy, a spectral flatness measure, power-law fits, and the number of
% crossings of the spectrum at various amplitude thresholds.
%
% Note there are two distinct senses of 'moment' here, which the field
% names can obscure: mom3 etc. are moments of the distribution of power
% *values* (order-agnostic in frequency -- shuffling the spectrum
% bin-for-bin leaves them unchanged), whereas specSkew/specKurt are
% moments of *frequency itself*, weighted by power (so they describe
% where in frequency the power sits). See the note on the
% power-weighted-moment block below.
%
% A linear-scale power spectrum is typically extremely heavy-tailed (cf.
% the peak-detection notes below), so several statistics that compare
% values or thresholds across the whole spectrum -- spectral
% autocorrelation (ac1/ac2/tau), quantiles (q25/median/q75), skewness
% (mom3), and the per-band stationarity measures (statav2_s/statav5_s)
% -- end up dominated by wherever the single largest value sits, rather
% than reflecting the spectrum's broader shape. Each of these has a
% log-domain companion (logac1/logac2/logtau, logq25/logmedian/logq75,
% logmom3, logstatav2_s/logstatav5_s) that applies the same statistic to
% logS instead: empirically validated to capture genuinely different
% information (|r| < 0.8 against the linear counterpart on two
% independent real datasets, mostly 0.2-0.5) rather than being a
% near-duplicate.
%
% Two families originally had both a linear and log-domain version --
% crossing counts (ncross) and the per-band mean-stationarity measures
% (statav2_m/statav5_m) -- but were cut back to log-only after a
% broader check: correlating every candidate against all ~7300 other
% hctsa features (not just a hand-picked "trivial" set) on 1000 real
% series showed the *linear* versions consistently had the lowest
% correlation with anything else in hctsa (R^2 ~ 0.1-0.3, vs ~0.2-0.6 for
% their log counterparts, which tended to land near existing
% autocorrelation/model-fit features -- i.e. real structure). Low
% redundancy with existing features is necessary but not sufficient for
% a feature to be useful (it's equally consistent with "just noise"), so
% this was a judgment call, not a statistical proof -- made here in
% favor of dropping the linear versions, consistent with them not having
% been observed to be useful in practice. ncross_f*/statav2_m/statav5_m
% are therefore no longer computed at all; only their logS counterparts
% remain.

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
%% Check that a Curve-Fitting Toolbox license is available:
% ------------------------------------------------------------------------------
BF_CheckToolbox('curve_fitting_toolbox');

% ------------------------------------------------------------------------------
% Check inputs, set defaults:
% ------------------------------------------------------------------------------
if size(y, 2) > size(y, 1)
	y = y'; % Time series must be a column vector
end
if nargin < 2 || isempty(psdMeth)
	psdMeth = 'fft'; % fft by default
end
if nargin < 3 || isempty(windowType)
	windowType = 'hamming'; % Hamming window by default
end
if nargin < 4
	nf = [];
end
if nargin < 5 || isempty(doLogAbs)
	doLogAbs = false;
end

if doLogAbs % a boolean
	% Analyze the spectrum of logarithmic absolute deviations
	y = log(abs(y));
end

doPlot = false; % plot outputs
Ny = length(y); % time-series length

% -------------------------------------------------------------------------------
% Set window (for periodogram and welch):
% -------------------------------------------------------------------------------
if ismember(psdMeth, {'periodogram', 'welch'})
	% Welch's method needs a window *shorter* than the full series so pwelch
	% actually splits it into multiple overlapping segments to average over --
	% that segment-averaging (trading frequency resolution for reduced
	% estimate variance) is the entire point of the method. A window spanning
	% the whole series (as used below for periodogram, a genuine
	% single-segment, whole-series estimate) would silently collapse Welch's
	% method to a single "segment" with no averaging at all, making it
	% numerically near-indistinguishable from an unwindowed periodogram: an
	% actual bug found and fixed here (validated on HCTSA_Empirical1000.mat,
	% where the un-fixed version had r=1.000 with the 'fft' method's output
	% across 97.6% of fields). noverlap is left at pwelch's own default (50%
	% of window length) below.
	if strcmp(psdMeth, 'welch')
		% winLength used to be Ny/4, which keeps the segment COUNT fixed
		% (~4, ~7 with 50% overlap) but lets the segment LENGTH -- and with
		% it the intrinsic frequency resolution and smoothing bandwidth --
		% grow without bound as the series lengthens. That defeats the
		% actual convergence guarantee of Welch's method (per-bin variance
		% should shrink as more segments are averaged over a *fixed*
		% resolution) and was the direct cause of several fields drifting
		% with N: verified with an AR(2) resonance of known -3dB width
		% (0.42 rad/sample) that the old scheme's measured maxWidth decayed
		% monotonically and without bound (0.23 -> 0.008 rad/sample, N =
		% 200 -> 6400, a 30x error at the largest length, actively getting
		% worse with more data) while capping winLength at a fixed absolute
		% length converges instead (0.15-0.20 across the same range).
		% 256 matches both scipy.signal.welch's and pwelch's own default
		% segment length (see the nfft note below), so above Ny ~ 1024 this
		% now behaves like a textbook Welch estimate: resolution fixed,
		% more segments averaged as N grows, variance shrinks with N.
		% Below that it is unchanged from before (round(Ny/4) is smaller
		% than 256), so short-series values are untouched.
		winLength = max(min(256, round(Ny / 4)), 16);
	else
		winLength = Ny;
	end
	switch windowType % method to use for the window
		case 'none'
			window = [];
		case 'hamming'
			window = hamming(winLength);
		case 'hann'
			window = hann(winLength);
		case 'bartlett'
			window = bartlett(winLength);
		case 'boxcar'
			window = boxcar(winLength);
		case 'rect'
			window = rectwin(winLength);
		otherwise
			% There are other options, but these aren't implemented here
			error('Unknown window, ''%s''', windowType);
	end
end

% ------------------------------------------------------------------------------
% Compute the Fourier Transform
% ------------------------------------------------------------------------------
switch psdMeth
	case 'periodogram'
		if isempty(nf)
			% (2) Estimate the spectrum
			[S, w] = periodogram(y, window);
		else
			w = linspace(0, pi, nf);
			[S, w] = periodogram(y, window, w);
		end

	case 'fft'
		% Fast Fourier Transform
		Fs = 1; % sampling frequency
		NFFT = 2^nextpow2(Ny);
		f = Fs / 2 * linspace(0, 1, NFFT / 2 + 1); % frequency
		w = 2 * pi * f'; % angular frequency (as column vector)
		S = fft(y, NFFT); % Fourier Transform
		S = 2 * abs(S(1:NFFT / 2 + 1)).^2 / Ny; % single-sided power spectral density
		S = S / (2 * pi); % convert to angular frequency space

	case 'welch'
		% Welch power spectral density estimate:
		Fs = 1; % sampling frequency
		% nfft used to be tied to Ny (2^nextpow2(Ny), same as the 'fft'
		% branch above), independently of winLength -- so even after
		% capping winLength, the frequency GRID kept getting finer with N
		% while the actual smoothing bandwidth (set by winLength) stayed
		% fixed. That mismatch showed up as its own N-dependence in
		% bin-lag statistics (ac1/ac2/logac1/logac2, which correlate
		% adjacent *bins*, not a physical frequency lag) and inflated
		% several bin-count sums (spect_shann_ent, the poly/saturating
		% cumsum fits) that are supposed to describe the spectrum's shape,
		% not its number of samples -- verified: e.g. fpolysat_a's eta^2
		% against length dropped from 0.98 to 0.06 once nfft stopped
		% growing with Ny. Passing [] lets pwelch fall back to its own
		% documented default (max(256, 2^nextpow2(winLength)) --
		% confirmed empirically), i.e. resolution set by the window, not
		% the series -- the textbook convention this file's manual
		% override had been silently bypassing.
		[S, f] = pwelch(y, window, [], [], Fs);
		w = 2 * pi * f'; % angular frequency
		S = S / (2 * pi); % adjust so that area remains normalized in angular frequency space

	otherwise
		error('Unknown spectral estimation method ''%s''', psdMeth);
end

if ~any(isfinite(S)) % no finite values in the power spectrum
	% This time series must be really weird -- return NaN (unsuitable operation)...
	warning('NaN in power spectrum? A weird time series.');
	out = NaN; return
end

% Ensure both w and S are row vectors:
if size(S, 1) > size(S, 2)
	S = S';
end
if size(w, 1) > size(w, 2)
	w = w';
end

if doPlot
	figure('color', 'w')
	plot(w, S, '.-k'); % plot the spectrum
	% Area under S should sum to 1 if a power spectral density estimate:
	title(sprintf('Area under psd curve = %.1f (= %.1f)', sum(S * (w(2) - w(1))), var(y)));
end

N = length(S); % = length(w)
logS = log(S);
dw = w(2) - w(1); % spacing increment in w

% ------------------------------------------------------------------------------
% Simple measures of the power spectrum
% ------------------------------------------------------------------------------

% -------------------------------------------------------------------------------
% Peaks:
% -------------------------------------------------------------------------------
% Maximum, and max peak width:
[out.maxS, i_maxS] = max(S);
out.maxw = w(i_maxS);

% Half-power (-3 dB) bandwidth of the dominant peak: the frequency interval
% around the maximum over which the spectrum stays above half its peak value.
%
% This previously thresholded against out.maxS itself. Since maxS *is* the
% maximum, every neighbouring bin satisfies S < maxS, so the search always
% terminated on the immediately adjacent bins and the field was identically
% 2*dw = 4*pi/nfft -- a pure function of the transform length carrying no
% information about the data (verified against an AR(2) resonance of known
% intrinsic width, where it still returned exactly 2*dw at every N).
halfPower = out.maxS / 2;
iUpper = i_maxS + find(S(i_maxS + 1:end) < halfPower, 1, 'first');
if isempty(iUpper) % never drops below half power above the peak
	iUpper = length(S);
end
iLower = find(S(1:i_maxS - 1) < halfPower, 1, 'last');
if isempty(iLower) % never drops below half power below the peak
	iLower = 1;
end
out.maxWidth = w(iUpper) - w(iLower);

% Characterize all peaks using findpeaks function, run on logS rather than S:
% a linear-scale power spectrum is typically extremely heavy-tailed (a single
% dominant peak can be >100x the mean level -- cf. the SP_Summaries curation
% notes), so prominence-based peak detection on raw S systematically buries
% smaller-but-genuine peaks (e.g. harmonics) under the dominant one, and
% findpeaks' fixed prominence thresholds below are only meaningfully
% calibrated on the log scale. This also matches Matlab's own convention:
% periodogram/pwelch plot in dB (log power) by default, not linear power.
%
% These prominence thresholds are calibrated specifically for a
% Welch-smoothed logS (i.e. psdMeth = 'welch'): a raw single-segment
% periodogram/fft estimate is a statistically inconsistent spectral
% estimator (per-bin variance doesn't shrink with more data), so its
% logS has a much noisier, less-calibrated peak structure -- validated via
% a white-noise null model (Welch-smoothed: 99th-percentile max spurious
% log-prominence ~3.5 across series lengths 500-10000; unsmoothed: ~12,
% overlapping real harmonic peaks almost entirely). For 'fft'/'periodogram',
% these fields are therefore left uncalibrated (computed for code-path
% uniformity only) and are deliberately NOT registered as hctsa features
% for those two psdMeth variants -- only for 'welch'.
minDist_w = 0.02;
ptsPerw = length(S) / pi;
minPkDist = ceil(minDist_w * ptsPerw);
[pkHeight, pkLoc, pkWidth, pkProm] = findpeaks(logS, 'SortStr', 'descend', 'minPeakDistance', minPkDist);
pkWidth = pkWidth / ptsPerw;
pkLoc = pkLoc / ptsPerw;

% Characterize mean peak prominence (thresholds in log-power units; see note above)
%
% numPeaks and numPeaks_overmean are computed (numPeaks is needed below) but are
% deliberately NOT registered as hctsa features, for the same reason the
% fft/periodogram peak fields above are not: they apply no prominence threshold,
% so they count noise-floor ripple and therefore measure spectral resolution
% rather than the data. On white noise numPeaks runs 16 -> 105 across
% N = 200..6400, and numPeaks_overmean returns essentially identical values for
% white noise (7 -> 53) and for a clear 3-tone signal (6 -> 43) -- i.e. no
% discriminative content. Normalizing by the number of resolvable peak slots
% does not fix this (0.12 -> 0.69 over the same range). The prominence-
% thresholded counts below are the calibrated, length-stable versions:
% numPromPeaks_3 returns exactly 0 for white noise and exactly 3 for the 3-tone
% signal at every length tested.
numPeaks = length(pkHeight); % local only: needed for the peakPower_* fields below, not itself an output
out.numPromPeaks_3 = sum(pkProm > 3); % number of peaks with log-prominence of at least 3 (~90-95th pctile of the white-noise null)
out.numPromPeaks_5 = sum(pkProm > 5); % number of peaks with log-prominence of at least 5 (clearly above the null's observed max of ~3.8)
out.numPromPeaks_8 = sum(pkProm > 8); % number of peaks with log-prominence of at least 8 (matches the weakest peak in a validated harmonic-series test signal)
out.meanProm_5 = mean(pkProm(pkProm > 5)); % mean peak prominence of those with log-prominence of at least 5

out.meanPeakWidth_prom5 = mean(pkWidth(pkProm > 5)); % mean peak width of peaks with log-prominence of at least 5
out.width_weighted_prom = sum(pkWidth .* pkProm) / sum(pkProm);

% Power in top N peaks:
nn = @(x) 1:min(x, numPeaks);
out.peakPower_2 = sum(pkHeight(nn(2)) .* pkWidth(nn(2)));
out.peakPower_5 = sum(pkHeight(nn(5)) .* pkWidth(nn(5)));
out.peakPower_prom5 = sum(pkHeight(pkProm > 5) .* pkWidth(pkProm > 5)); % power in peaks with log-prominence of at least 5
out.w_weighted_peak_prom = sum(pkLoc .* pkProm) / sum(pkProm); % where are prominent peaks located on average (weighted by prominence)

% Number of peaks required to get to 50% of power in peaks
peakPower = pkHeight .* pkWidth;
if isempty(peakPower) % no peaks at all (e.g. a monotonic power spectrum)
    out.numPeaks_50power = NaN;
    out.peakpower_1 = NaN;
else
    out.numPeaks_50power = find(cumsum(peakPower) > 0.5 * sum(peakPower), 1, 'first');
    out.peakpower_1 = peakPower(1) / sum(peakPower);
end

% -------------------------------------------------------------------------------
% Distribution
% -------------------------------------------------------------------------------
% Quantiles:
out.iqr = iqr(S);
out.logiqr = iqr(logS);
out.q25 = quantile(S, 0.25);
out.median = median(S);
out.q75 = quantile(S, 0.75);
% Log-domain quantiles: q25/median/q75 above only have a linear version,
% even though the iqr computed from them already has a logiqr companion
% -- add the missing log-domain counterparts for consistency.
out.logq25 = quantile(logS, 0.25);
out.logmedian = median(logS);
out.logq75 = quantile(logS, 0.75);

% Moments:
out.std = std(S);
out.stdlog = log(out.std);
out.logstd = std(logS);
out.mom3 = DN_Moments(S, 3, true);
out.logmom3 = DN_Moments(logS, 3, true);

% Autocorrelation of amplitude spectrum:
autoCorrs_S = CO_AutoCorr(S, 1:4, 'Fourier');
out.ac1 = autoCorrs_S(1);
out.ac2 = autoCorrs_S(2);
% CO_FirstCrossing returns a lag measured in *bins*, not physical
% frequency -- and the number of bins (NFFT) scales with the input
% series length N, independent of any real spectral structure. Left
% unconverted, tau is partly just re-deriving N (checked empirically:
% R^2 = 0.70 against hctsa's own `length` feature across 1000 real
% series, by far the single biggest external correlate found). Multiply
% by dw (the bin spacing in angular frequency) to convert to a genuine,
% resolution-independent physical unit: doubling NFFT halves dw but
% doubles the bin-count of any fixed physical correlation width, so the
% product is invariant to N, leaving only the spectrum's actual
% correlation structure.
out.tau = CO_FirstCrossing(S, 'ac', 0, 'continuous') * dw;

% Same, on logS: Pearson autocorrelation on the raw (heavy-tailed) linear
% spectrum is effectively dominated by however far the single largest
% value sits from the rest -- one point can carry most of the
% correlation's covariance term. logS compresses that dominance so the
% autocorrelation instead reflects genuine bin-to-bin shape structure.
autoCorrs_logS = CO_AutoCorr(logS, 1:4, 'Fourier');
out.logac1 = autoCorrs_logS(1);
out.logac2 = autoCorrs_logS(2);
out.logtau = CO_FirstCrossing(logS, 'ac', 0, 'continuous') * dw; % see tau's N-independence note above

% ------------------------------------------------------------------------------
% Shape of cumulative sum curve
% ------------------------------------------------------------------------------
csS = cumsum(S);

f_frac_w_max = @(f) w(find(csS >= csS(end) * f, 1, 'first'));

% At what frequency is csS a fraction p of its maximum?
out.wmax_5 = f_frac_w_max(0.05);
out.wmax_10 = f_frac_w_max(0.1);
out.wmax_25 = f_frac_w_max(0.25);
out.centroid = f_frac_w_max(0.5);
out.wmax_75 = f_frac_w_max(0.75);
out.wmax_90 = f_frac_w_max(0.9);
out.wmax_95 = f_frac_w_max(0.95);
out.wmax_99 = f_frac_w_max(0.99);

% ------------------------------------------------------------------------------
% Power-weighted moments of the frequency distribution
% ------------------------------------------------------------------------------
% Note these are a genuinely different family from both of the
% frequency-summarizing families above and the moments of S computed
% earlier, in a way the naming can obscure:
%   - `centroid` above is the *median* frequency (where csS reaches 50%),
%     NOT the power-weighted mean frequency computed here. It keeps that
%     (strictly speaking, misnamed) name because it is one of the
%     canonical catch22 features (SP_Summaries_welch_rect_centroid),
%     externally mirrored in the catch22 C/Python implementations -- so
%     the name carries a compatibility contract and must not change.
%   - `mom3` (skewness of S) summarizes the distribution of power
%     *values*, which is entirely order-agnostic in frequency: shuffle
%     the spectrum bin-for-bin and mom3 is unchanged. It says nothing
%     about *where* in frequency the power sits.
% The moments below instead treat the (non-negative) spectrum as a
% weighting over frequency and take moments of frequency itself -- the
% standard spectral centroid/spread/skewness/kurtosis (cf. the MPEG-7
% audio descriptors). Relative to the wmax_* quantiles, which describe
% the same frequency distribution non-parametrically, these weight the
% tails much more heavily: a small amount of power far out in frequency
% moves specSkew/specKurt substantially and wmax_90 not at all.
%
% There is deliberately no log-domain companion here (unlike most of the
% families above): S enters as a *weight*, not as a value being
% summarized, and logS is negative wherever S < 1, so log-weighting
% would produce meaningless (signed) 'weights'.
%
% On redundancy: checked on 1000 real series, the two low-order moments
% are substantially re-descriptions of the wmax_* quantiles already
% computed above -- specCentroid vs wmax_75 R^2 = 0.91 (and vs the
% median-frequency `centroid`, 0.88), specSpread vs the 10-90% width
% w10_90 R^2 = 0.92. That is unsurprising (a mean tracks a median, a
% standard deviation tracks an inter-quantile width). They are kept and
% registered anyway, deliberately: specCentroid and specSpread are the
% textbook spectral centroid and bandwidth (the MPEG-7 audio
% descriptors), so they are far easier to interpret and to motivate
% theoretically than the quantile-width proxies that happen to correlate
% with them -- worth a little redundancy. specSkew and specKurt are
% additionally non-redundant on the numbers: worst-case R^2 of 0.20 and
% 0.08 against any of those quantiles, since they respond to
% far-out-in-frequency tail power the quantiles are by construction
% insensitive to.
Spos = max(S, 0); % guard against any tiny negative values from the estimator
sumS = sum(Spos);
if sumS > 0
    pw = Spos / sumS; % normalized weighting over frequency
    out.specCentroid = sum(pw .* w);
    wDev = w - out.specCentroid;
    specVar = sum(pw .* wDev.^2);
    out.specSpread = sqrt(specVar);
    if specVar > 0
        out.specSkew = sum(pw .* wDev.^3) / specVar^(3/2);
        out.specKurt = sum(pw .* wDev.^4) / specVar^2;
    else % all power in a single frequency bin: shape undefined
        out.specSkew = NaN;
        out.specKurt = NaN;
    end
else % no power anywhere
    out.specCentroid = NaN;
    out.specSpread = NaN;
    out.specSkew = NaN;
    out.specKurt = NaN;
end

% ------------------------------------------------------------------------------
% Fit some functions to this cumulative sum:
% ------------------------------------------------------------------------------
% (i) Quadratic
[c, gof] = fit(w', csS', 'poly2');
out.fpoly2csS_p1 = c.p1;
out.fpoly2csS_p2 = c.p2;
out.fpoly2csS_p3 = c.p3;
out.fpoly2_sse = gof.sse;
out.fpoly2_r2 = gof.rsquare;
out.fpoly2_rmse = gof.rmse;

% (ii) Fit polysat a*x^2/(b+x^2) (has zero derivative at zero, though)
s = fitoptions('Method', 'NonlinearLeastSquares', 'StartPoint', [csS(end), 100]);
f = fittype('a*x^2/(b+x^2)', 'independent', 'x', 'options', s); % set 'a' from maximum
[c, gof] = fit(w', csS', f);
out.fpolysat_a = c.a;
out.fpolysat_b = c.b; % this is important
out.fpolysat_r2 = gof.rsquare; % this is more important!
out.fpolysat_rmse = gof.rmse;

% ------------------------------------------------------------------------------
% Shannon spectral entropy
% ------------------------------------------------------------------------------
Hshann = -S .* log(S); % Shannon function
out.spect_shann_ent = sum(Hshann);
out.spect_shann_ent_norm = mean(Hshann);

% ------------------------------------------------------------------------------
% "Spectral Flatness Measure"
% ------------------------------------------------------------------------------
% which is given in dB as 10 log_10(gm/am) where gm is the geometric mean and am
% is the arithmetic mean of the power spectral density
out.sfm = 10 * log10(geomean(S) / mean(S));

% ------------------------------------------------------------------------------
% Areas under power spectrum
% ------------------------------------------------------------------------------
% Area up to peak: (may be more appropriate in squared log units?)
out.areatopeak = sum(S(1:i_maxS)) * dw;
out.ylogareatopeak = sum(logS(1:i_maxS)) * dw; % (semilogy)
% out.logareatopeak=sum(logS(1:i_maxS).*diff(logw(1:i_maxS+1)));

% ------------------------------------------------------------------------------
%% Robust (e.g., iteratively re-weighted least squares) linear fits to log-log
%   plots
% ------------------------------------------------------------------------------

% Use the local function, giveMeRobustStats, which adds robust
% stats to the out structure.

% Suppress rank deficient warnings for this section:
warning('off', 'stats:robustfit:RankDeficient')

% (1): Across full range
r_all = (w > 0); % avoid -Inf for log(0) when w = 0;
out = giveMeRobustStats(log(w(r_all)), log(S(r_all)), 'linfitloglog_all', out, {'a1','a2','sigrat','sigma','sea1'});

% (2): First half (low frequency)
r_lf = (w > 0); % w(1) = 0 -> log(0) = -Inf
r_lf(floor(N / 2) + 1:end) = 0; % remove second half of angular frequenciesf
out = giveMeRobustStats(log(w(r_lf)), log(S(r_lf)), 'linfitloglog_lf', out, {'a2'});

% (3): Second half (high frequency)
r_hf = floor(N / 2) + 1:N;
out = giveMeRobustStats(log(w(r_hf)), log(S(r_hf)), 'linfitloglog_hf', out, {'a1','a2','sigrat','sigma','sea1'});

% (4): Middle half (mid-frequencies)
r_mf = round(N / 4):round(N * 3 / 4);
out = giveMeRobustStats(log(w(r_mf)), log(S(r_mf)), 'linfitloglog_mf', out, {'a2'});

% (5) Fit linear to semilog plot (across full range)
out = giveMeRobustStats(w, log(S), 'linfitsemilog_all', out, {'a1','sigrat','sigma','sea1'});

% Turn the rank-deficient warnings back on
warning('on', 'stats:robustfit:RankDeficient')

% ------------------------------------------------------------------------------
%% Power in specific frequency bands
% ------------------------------------------------------------------------------
% *** DO THIS BY BUFFER COMMAND: AND WHILE AT IT LOOK AT STATIONARITY OF
% POWER SPECTRUM....

% 2 bands
split = buffer(S, floor(N / 2));
if size(split, 2) > 2, split = split(:, 1:2); end
out.area_2_1 = sum(split(:, 1)) * dw;
out.logarea_2_1 = sum(log(split(:, 1))) * dw;
out.area_2_2 = sum(split(:, 2)) * dw;
out.logarea_2_2 = sum(log(split(:, 2))) * dw;
out.statav2_s = std(std(split)) / std(S);
% Same on logS: whichever band happens to contain the dominant peak
% swamps std(mean(split))/std(std(split)) on the raw linear spectrum, so
% this is essentially "which band has the peak" rather than a nuanced
% measure of stationarity across the spectrum's shape. (The mean-based
% version, statav2_m, was dropped entirely -- see the note at the top of
% this file; std(mean(split))/std(S) checked out as the least-redundant-
% with-anything candidate in the same cross-hctsa audit that flagged
% ncross, i.e. noise rather than novel structure.)
splitLog = buffer(logS, floor(N / 2));
if size(splitLog, 2) > 2, splitLog = splitLog(:, 1:2); end
out.logstatav2_m = std(mean(splitLog)) / std(logS);
out.logstatav2_s = std(std(splitLog)) / std(logS);

% 5 bands
split = buffer(S, floor(N / 5));
if size(split, 2) > 5, split = split(:, 1:5); end
out.area_5_1 = sum(split(:, 1)) * dw;
out.logarea_5_1 = sum(log(split(:, 1))) * dw;
out.area_5_2 = sum(split(:, 2)) * dw;
out.logarea_5_2 = sum(log(split(:, 2))) * dw;
out.area_5_3 = sum(split(:, 3)) * dw;
out.logarea_5_3 = sum(log(split(:, 3))) * dw;
out.area_5_4 = sum(split(:, 4)) * dw;
out.logarea_5_4 = sum(log(split(:, 4))) * dw;
out.area_5_5 = sum(split(:, 5)) * dw;
out.logarea_5_5 = sum(log(split(:, 5))) * dw;
out.statav5_s = std(std(split)) / std(S);
splitLog = buffer(logS, floor(N / 5));
if size(splitLog, 2) > 5, splitLog = splitLog(:, 1:5); end
out.logstatav5_m = std(mean(splitLog)) / std(logS);
out.logstatav5_s = std(std(splitLog)) / std(logS);

% ------------------------------------------------------------------------------
% Count crossings:
% Get a horizontal line and count the number of crossings with the power spectrum
% ------------------------------------------------------------------------------
% A threshold that's a fixed fraction of max(S) sits far above the noise
% floor whenever there's one dominant peak (the same heavy-tailed-
% spectrum problem that motivated moving peak detection to logS), so
% crossing counts on the raw linear spectrum end up governed almost
% entirely by the shape immediately around that one peak -- and,
% checked empirically against ~7300 other hctsa features on 1000 real
% series, ended up the least redundant with anything else in hctsa of
% any candidate checked (R^2 ~ 0.1-0.3), consistent with that being
% noise rather than novel structure. Dropped in favor of the logS-only
% version below, where differences correspond to power *ratios* (dB), so
% a threshold set as a fraction of the full [min(logS), max(logS)] range
% is a genuinely relative "how far up from the noise floor" level, not
% just "how close to the single biggest value":
logRange = max(logS) - min(logS);
ncrossfn_rel_log = @(f) sum(BF_SignChange(logS - (min(logS) + f * logRange)));

out.ncross_log_f05 = ncrossfn_rel_log(0.05);
out.ncross_log_f10 = ncrossfn_rel_log(0.1);
out.ncross_log_f20 = ncrossfn_rel_log(0.2);
out.ncross_log_f50 = ncrossfn_rel_log(0.5);

% -------------------------------------------------------------------------------
% function mel = w2mel(w) % convert to mel spectrum
%     mel = 1127*log(w/(1400*pi)+1);
% end

function out = giveMeRobustStats(xData, yData, textID, out, whichStats)
	% Add statistics to the output structure from a robust linear fit
	% between xData and yData.
	%
	% whichStats selects which of the six available statistics to emit.
	% It exists because the full set was previously emitted for every fit
	% range, but only a subset of those were ever registered as hctsa
	% features -- for the 'lf' and 'mf' ranges only the gradient was --
	% leaving 14 fields computed and returned but never used by anything.

	% Perform the fit:
	[a, stats] = robustfit(xData, yData);

	% Add the requested statistics to the output structure:
	if ismember('a1', whichStats)
		out.(sprintf('%s_a1', textID)) = a(1); % robust intercept
	end
	if ismember('a2', whichStats)
		out.(sprintf('%s_a2', textID)) = a(2); % robust gradient
	end
	if ismember('sigrat', whichStats)
		% ratio of sigma estimates between ordinary least squares (ols) and the robust fit:
		out.(sprintf('%s_sigrat', textID)) = stats.ols_s / stats.robust_s;
	end
	if ismember('sigma', whichStats)
		% esimate of sigma as the larger of robust_s and a weighted average of ols_s and robust_s:
		out.(sprintf('%s_sigma', textID)) = stats.s;
	end
	if ismember('sea1', whichStats)
		out.(sprintf('%s_sea1', textID)) = stats.se(1); % standard error of 1st coefficient estimate
	end
	if ismember('sea2', whichStats)
		out.(sprintf('%s_sea2', textID)) = stats.se(2); % standard error of 2nd coefficient estimate
	end
end

end
