function out = PP_Compare(y, detrndmeth)
% PP_Compare   How time-series properties change after a preprocessing step.
%
% Applies a given preprocessing transformation (detrending, differencing,
% filtering or resampling) to the time series, z-scores the original and the
% processed series, and returns the change in each of a set of statistics from its
% value for the original to its value for the processed series. The statistics compare
% stationarity measures (StatAv, and the variation of the local mean and of the
% local standard deviation across windows), distributional fits (a Gaussian fit
% to the kernel-smoothed distribution, and the discrepancy from a fitted normal
% distribution) and the effect of trimming outliers.
%
% The change is the difference (processed minus original) for statistics that can be
% negative or zero, and the normalized difference (processed - original) /
% (processed + original) for positive statistics. The latter is between -1 and 1, is 0
% when nothing changes, and equals tanh(log(processed / original) / 2), a bounded
% version of the log-ratio that stays finite when the original value is near 0, where
% a plain ratio is unstable.
%
% ---INPUTS:
% y, the input time series
% detrndmeth, the preprocessing to apply (default: 'medianf3'):
%       'poly<n>': remove a polynomial of order n = 1-9 (Curve Fitting Toolbox),
%           e.g., 'poly1', a linear detrending
%       'sin<n>': remove a constant plus a sum of n = 1-8 sinusoids
%           a1*sin(2*pi*f1*t+c1) + ..., fitted jointly by least squares, with frequencies
%           searched deterministically between 1/N and 1/2 - 1/(2N) cycles per sample
%           (BF_FitSinusoids), e.g., 'sin1': this removes oscillations, and a component
%           with less than one full cycle in the series is a trend, not an oscillation
%           (trends are removed by 'poly' and 'spline')
%       'spline<npieces><order>': remove a least-squares spline fitted with
%           spap2 (Spline Toolbox) with the given number of polynomial pieces and
%           spline order, e.g., 'spline24', a cubic spline with 2 pieces
%       'diff<n>': take n successive differences, e.g., 'diff3'
%       'medianf<n>': a running median filter of length n (centered; the end values
%           are repeated beyond the ends), e.g., 'medianf3'
%       'rav<n>': a running mean filter of length n (filter), e.g., 'rav5'
%       'resample_<p>_<q>': resample the series by the ratio p/q (resample); e.g.,
%           'resample_1_2' halves the length and 'resample_10_1' multiplies it by 10
%       'logr': log returns (positive data only; otherwise NaN)
%       'boxcox': a Box-Cox transformation (positive data only; otherwise NaN)
%
% ---OUTPUTS: the change, from the original to the processed series, of each of these
% statistics (all of the series are z-scored first):
% Normalized differences (positive statistics):
% statav2, StatAv with 2 segments (SY_StatAv)
% swms2_2, swms5_1, swms10_1, the standard deviation of the window means across
%       windows (SY_SlidingWindow 'mean'), with 2 windows overlapping by half, and
%       5 and 10 non-overlapping windows
% swss2_1, swss5_1, swss10_1, the same for the window standard deviations
%       (SY_SlidingWindow 'std')
% kscn_olapint, the overlap integral of the kernel-smoothed distribution with the
%       best-fitting normal (DN_CompareKSFit)
% olbt_s5, the standard deviation after trimming the 5% most extreme values at each
%       end, relative to that of the full series (DN_OutlierTest)
% Differences (statistics that can be negative or zero):
% gauss1_kd_r2, gauss1_kd_resAC1, gauss1_kd_resrunsz, the R^2, the lag-1
%       autocorrelation of the residuals, and the runs-test z-statistic of the residuals
%       of a Gaussian fit to the kernel-smoothed distribution (DN_SimpleFit)
% kscn_peaksepy, kscn_peaksepx, kscn_relent, the peak separation in height and in
%       position, and the relative entropy, of the kernel-smoothed distribution
%       against the best-fitting normal (DN_CompareKSFit)
% olbt_m2, olbt_m5, the mean after trimming the 2% and 5% most extreme values at
%       each end (DN_OutlierTest)
% A scalar NaN is returned if the processed series has no variation beyond
% floating-point round-off, as statistics of round-off are not meaningful: this is the
% case for the residual of an exact fit (e.g., a polynomial or a sum of sinusoids that
% the series is), or the differences of an exact linear trend. The criterion is a
% standard deviation of at most 10*N*eps times the root-mean-square of the series the
% processing acts on (its length N; the log of the series for 'logr', and the
% transformed series itself for 'boxcox'), a bound on the rounding error of the
% processing (sums of up to N terms, and the phase of sinusoids evaluated at times up
% to N). Small but genuine variation, such as the quantization noise left by an exact
% fit to 24-bit data (about 1e-7 of the amplitude), is computed. NaN is also returned
% for 'logr' and 'boxcox' if the data are not all positive.

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
%% Check inputs, set default:
if nargin < 2 || isempty(detrndmeth)
	detrndmeth = 'medianf3'; % median filter by default
end

% -------------------------------------------------------------------------------
% Preparation
N = length(y); % time-series length
r = (1:N)'; % the time-range over which to fit
yRef = y; % the series the processing acts on (its magnitude sets the round-off below)

% -------------------------------------------------------------------------------
%% APPLY PREPROCESSINGS
% ------------------------------------------------------------------------------
% DETRENDINGS:
% Do the detrending; converting from y (raw) to y_d (detrended) by
% subtracting some fit y_fit

% 1) Polynomial detrend
% starts with 'poly' and ends with integer from 1--9
if length(detrndmeth) == 5 && strcmp(detrndmeth(1:4), 'poly') && ~isnan(str2double(detrndmeth(5)))

	% Check a curve-fitting toolbox license is available:
	BF_CheckToolbox('curve_fitting_toolbox');

	% (time centered and scaled, so that high orders are well conditioned)
	[cfun, ~] = fit(r, y, detrndmeth, 'Normalize', 'on');
	y_fit = feval(cfun, r);
	y_d = y - y_fit;

	% 2) Seasonal detrend
elseif length(detrndmeth) == 4 && strcmp(detrndmeth(1:3), 'sin') && ~isnan(str2double(detrndmeth(4))) && ~strcmp(detrndmeth(4), '9')

	% (a constant is fitted jointly, so that a constant plus sinusoids is removed exactly:
	% the frequencies are bounded away from zero, so the sinusoids cannot absorb a mean)
	numSin = str2double(detrndmeth(4));
	if N <= 3*numSin + 1
		out = NaN; % too short to fit this many sinusoids and a constant
		return
	end
	y_d = y - BF_FitSinusoids(y, numSin, true);

	% 3) Spline detrend
elseif length(detrndmeth) == 8 && strcmp(detrndmeth(1:6), 'spline') && ~isnan(str2double(detrndmeth(7))) && ~isnan(str2double(detrndmeth(8)))
	nknots = str2double(detrndmeth(7));
	intp = str2double(detrndmeth(8));

	%% Check that a Curve-Fitting Toolbox license is available:
	BF_CheckToolbox('curve_fitting_toolbox');

	spline = spap2(nknots, intp, r, y); % just a single middle knot with cubic interpolants
	y_spl = fnval(spline, 1:N); % evaluate at the 1:N time intervals
	y_d = y - y_spl';

	% 4) Differencing
elseif length(detrndmeth) == 5 && strcmp(detrndmeth(1:4), 'diff') && ~isnan(str2double(detrndmeth(5)))
	ndiffs = str2double(detrndmeth(5));
	y_d = diff(y, ndiffs); % difference the series n times

	% 5) Median Filter
elseif length(detrndmeth) > 7 && strcmp(detrndmeth(1:7), 'medianf') && ~isnan(str2double(detrndmeth(8:end)))
	n = str2double(detrndmeth(8:end)); % order of filtering
	% (the first and last values are repeated beyond the ends, so the filtered series is
	% not pulled toward zero there, which would depend on the offset of the data)
	lo = floor(n / 2); hi = n - 1 - lo;
	y_d = movmedian([repmat(y(1), lo, 1); y(:); repmat(y(end), hi, 1)], [lo, hi], 'Endpoints', 'discard');

	% 6) Running Average
elseif length(detrndmeth) > 3 && strcmp(detrndmeth(1:3), 'rav') && ~isnan(str2double(detrndmeth(4:end)))
	n = str2double(detrndmeth(4:end)); % the window size
	y_d = filter(ones(1, n) / n, 1, y);

	% 7) Resample
elseif length(detrndmeth) > 9 && strcmp(detrndmeth(1:9), 'resample_')
	% check a valid structure
	ss = textscan(detrndmeth, '%s%n%n', 'delimiter', '_');
	if ~isempty(ss{2})
		p = ss{2};
	else
		error('Invalid ''resample_p_q'' detrending specification: ''%s''', detrndmeth)
	end
	if ~isempty(ss{3})
		q = ss{3};
	else
		error('Invalid ''resample_p_q'' detrending specification: ''%s''', detrndmeth)
	end
	y_d = resample(y, p, q);

	% 8) Log Returns
elseif strcmp(detrndmeth, 'logr')
	if all(y > 0)
		yRef = log(y);
		y_d = diff(yRef);
	else
		out = NaN;
		return % return all NaNs
	end

	% 9) Box-Cox Transformation
elseif strcmp(detrndmeth, 'boxcox')
	% Requires a financial toolbox to run boxcox, check one is available:
	BF_CheckToolbox('financial_toolbox');

	if all(y > 0)
		y_d = boxcox(y);
		yRef = y_d; % a transformation, not a residual: only a constant result is degenerate
	else
		out = NaN;
		return % return all NaNs
	end
else
	error('Invalid detrending method ''%s''', detrndmeth)
end

% -------------------------------------------------------------------------------
%% Quick check that outputs are meaningful
% The processed series is round-off if its variation is within the rounding error of the
% processing, bounded relative to the magnitude (root-mean-square) of the values it acts
% on: a fitted or filtered value accumulates rounding over sums of up to N terms (at
% most about N*eps/2), and a sinusoid of frequency up to 1/2 evaluated at times up to N
% has a rounding error in its phase of up to about 1.5*pi*N*eps (three roundings of
% 2*pi*f*t), in both the data and the fit. Averaged over the series (the phase error
% grows with t), these are below about 8*N*eps of the root-mean-square, so 10*N*eps
% covers them.
if std(y_d) <= 10 * length(yRef) * eps * sqrt(mean(yRef.^2))
	out = NaN;
	return
end

% -------------------------------------------------------------------------------
% Statistical tests on original and processed time series
% z-score both (these metrics will need it, and not always done beforehand
% because of positive-only data, etc.)
y = zscore(y);
y_d = zscore(y_d);

% Changes from the original to the processed series: a difference for
% statistics that can be negative or near 0, and a normalized difference for
% positive statistics:
f_diff = @(proc, orig) proc - orig;
f_normDiff = @(proc, orig) SUB_normDiff(proc, orig);

% The statistics of the original series are the same for every preprocessing, so
% they are computed once per series (cached, as in BF_GetTau):
orig = SUB_origStats(y);

% 1) Stationarity

% (a) StatAv
out.statav2 = f_normDiff(SY_StatAv(y_d, 'seg', 2), orig.statav2);

% (b) Sliding window mean
out.swms2_2 = f_normDiff(SY_SlidingWindow(y_d, 'mean', 'std', 2, 2), orig.swms2_2);
out.swms5_1 = f_normDiff(SY_SlidingWindow(y_d, 'mean', 'std', 5, 1), orig.swms5_1);
out.swms10_1 = f_normDiff(SY_SlidingWindow(y_d, 'mean', 'std', 10, 1), orig.swms10_1);

% (c) Sliding window std
out.swss2_1 = f_normDiff(SY_SlidingWindow(y_d, 'std', 'std', 2, 1), orig.swss2_1);
out.swss5_1 = f_normDiff(SY_SlidingWindow(y_d, 'std', 'std', 5, 1), orig.swss5_1);
out.swss10_1 = f_normDiff(SY_SlidingWindow(y_d, 'std', 'std', 10, 1), orig.swss10_1);

% 2) Gaussianity
% (a) kernel density fit
me1 = DN_SimpleFit(y_d, 'gauss1', 0); % kernel density fit to 1-peak gaussian
me2 = orig.simpleFit;
if (~isstruct(me1) && isnan(me1)) || (~isstruct(me2) && isnan(me2))
	% fitting gaussian failed -- returns a NaN rather than a structure
	out.gauss1_kd_r2 = NaN;
	out.gauss1_kd_resAC1 = NaN;
	out.gauss1_kd_resrunsz = NaN;
else
	out.gauss1_kd_r2 = f_diff(me1.r2, me2.r2);
	out.gauss1_kd_resAC1 = f_diff(me1.resAC1, me2.resAC1);
	out.gauss1_kd_resrunsz = f_diff(me1.resrunsz, me2.resrunsz);
end

% (b) compare distribution to fitted normal distribution
me1 = DN_CompareKSFit(y_d, 'norm');
me2 = orig.compareKSFit;

out.kscn_peaksepy = f_diff(me1.peaksepy, me2.peaksepy);
out.kscn_peaksepx = f_diff(me1.peaksepx, me2.peaksepx);
out.kscn_olapint = f_normDiff(me1.olapint, me2.olapint);
out.kscn_relent = f_diff(me1.relent, me2.relent);

% 3) Outliers
out.olbt_m2 = f_diff(DN_OutlierTest(y_d, 2, 'mean'), orig.olbt_m2);
out.olbt_m5 = f_diff(DN_OutlierTest(y_d, 5, 'mean'), orig.olbt_m5);
out.olbt_s5 = f_normDiff(DN_OutlierTest(y_d, 5, 'std'), orig.olbt_s5);

% ------------------------------------------------------------------------------
function orig = SUB_origStats(y)
	% Statistics of the z-scored original series, shared by all preprocessings:
	% cache the last few results, keyed on the exact (z-scored) input series
	persistent cacheY cacheOrig
	maxCacheEntries = 4;
	if isempty(cacheY)
		cacheY = {};
		cacheOrig = {};
	end
	for ci = 1:numel(cacheY)
		if isequaln(y, cacheY{ci})
			orig = cacheOrig{ci};
			return
		end
	end
	orig.statav2 = SY_StatAv(y, 'seg', 2);
	orig.swms2_2 = SY_SlidingWindow(y, 'mean', 'std', 2, 2);
	orig.swms5_1 = SY_SlidingWindow(y, 'mean', 'std', 5, 1);
	orig.swms10_1 = SY_SlidingWindow(y, 'mean', 'std', 10, 1);
	orig.swss2_1 = SY_SlidingWindow(y, 'std', 'std', 2, 1);
	orig.swss5_1 = SY_SlidingWindow(y, 'std', 'std', 5, 1);
	orig.swss10_1 = SY_SlidingWindow(y, 'std', 'std', 10, 1);
	orig.simpleFit = DN_SimpleFit(y, 'gauss1', 0);
	orig.compareKSFit = DN_CompareKSFit(y, 'norm');
	orig.olbt_m2 = DN_OutlierTest(y, 2, 'mean');
	orig.olbt_m5 = DN_OutlierTest(y, 5, 'mean');
	orig.olbt_s5 = DN_OutlierTest(y, 5, 'std');
	cacheY = [{y}, cacheY];
	cacheOrig = [{orig}, cacheOrig];
	if numel(cacheY) > maxCacheEntries
		cacheY(maxCacheEntries+1:end) = [];
		cacheOrig(maxCacheEntries+1:end) = [];
	end
end

% ------------------------------------------------------------------------------
function nd = SUB_normDiff(proc, orig)
	% Normalized difference (proc - orig) / (proc + orig) of two positive numbers,
	% in [-1, 1]; 0 if both are 0
	if proc + orig == 0
		nd = 0;
	else
		nd = (proc - orig) / (proc + orig);
	end
end

end
