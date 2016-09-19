function out = PP_Compare(y,detrndmeth)
% PP_Compare    Compare how time-series properties change after pre-processing.
%
% Applies a given pre-processing transformation to the time series, and returns
% statistics on how various time-series properties change as a result.
%
% Inputs are structured in a clunky way, unfortunately:
%
%---INPUTS:
% y, the input time series
% detrndmeth, the method to use for detrending:
%      (i) 'poly': polynomial detrendings, both linear and quadratic. Can
%                  be of the following forms:
%            (a) polynomial of given order: 'poly1', 'poly2', 'poly3',
%                'poly4', 'poly5', 'poly6', 'poly7', 'poly8', 'poly9'
%            (b) fit best polynomial: 'polybest' determines 'best' by
%                            various tests (e.g., whiteness of residuals,
%                            etc.)
%            (c) 'fitstrong' only fits if a 'strong' trend.
%      (ii) 'sin': sinusoidal detrending with either one or two frequency
%                components,
%            (a) fit a sine series of a given order
%               Fits a form like: a1*sin(b1*x+c1) + a2*sin(b2*x+c2) + ...
%               Additional number determines how many terms to include in the
%               series: 'sin1', 'sin2', 'sin3', 'sin4', 'sin5', 'sin6', 'sin7',
%               'sin8'
%            (b) 'sin_st1': fit only if a strong trend (i.e., if the amplitudes
%                                           are above a given threshold)
%      (iii) 'spline': removes a least squares spline using Matlab's
%                      Spline Toolbox function spap2
%                      Input of the form 'spline<nknots><interpolant_order>'
%                      e.g., 'spline45' uses four knots and 5th order
%                      interpolants (Implements a least squares spline via the
%                      spline toolbox function spap2)
%      (iv) 'diff': takes incremental differences of the time series. Of form
%                 'diff<ndiff>' where ndiff is the number of differencings to
%                 perform. e.g., 'diff3' performs three recursive differences
%      (v) 'medianf': a running median filter using a given window lengths
%                   Of form 'medianf<n>' where n is the window length.
%                   e.g., 'medianf3' takes a running median using the median of
%                     every 3 consecutive values as a point in the filtered
%                     time series. Uses the Signal Processing Toolbox
%                     function medfilt1
%      (vi) 'rav': running mean filter of a given window length.
%                  Uses Matlab's filter function to perform a running
%                  average of order n. Of form 'rav<n>' where n is the order of
%                  the running average.
%      (vii) 'resample': resamples the data by a given ratio using the resample
%                        function in Matlab.
%                        Of form 'resample_<p>_<q>', where the ratio p/q is the
%                        new sampling rate e.g., 'resample_1_2' will downsample
%                        the signal by one half e.g., resample_10_1' will
%                        resample the signal to 10 times its original length
%      (viii) 'logr': takes log returns of the data. Only valid for positive
%                       data, else returns a NaN.
%      (ix) 'boxcox': makes a Box-Cox transformation of the data. Only valid for
%                     positive only data; otherwise returns a NaN.
%
% If multiple detrending methods are specified, they should be in a cell of
% strings; the methods will be executed in order, e.g., {'poly1','sin1'} does a
% linear polynomial then a simple one-frequency seasonal detrend (in that order)
%
%---OUTPUTS: include comparisons of stationarity and distributional measures
% between the original and transformed time series.

% ------------------------------------------------------------------------------
% Copyright (C) 2015, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013). DOI: 10.1098/rsif.2013.0048
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
%% Check inputs, set default
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(detrndmeth)
    detrndmeth = 'medianf'; % median filter by default
end

% ------------------------------------------------------------------------------
%% FOREPLAY
% ------------------------------------------------------------------------------
N = length(y); % time-series length
r = (1:N)'; % the t-range over which to fit

% ------------------------------------------------------------------------------
%% APPLY PREPROCESSINGS
% ------------------------------------------------------------------------------
% DETRENDINGS:
% Do the detrending; converting from y (raw) to y_d (detrended) by
% subtracting some fit y_fit

% 1) Polynomial detrend
% starts with 'poly' and ends with integer from 1--9
if length(detrndmeth) == 5 && strcmp(detrndmeth(1:4),'poly') && ~isempty(str2double(detrndmeth(5)))

    % Check a curve-fitting toolbox license is available:
    BF_CheckToolbox('curve_fitting_toolbox');

    [cfun, gof] = fit(r,y,detrndmeth);
    y_fit = feval(cfun,r);
    y_d = y - y_fit;

% 2) Seasonal detrend
elseif length(detrndmeth) == 4 && strcmp(detrndmeth(1:3),'sin') && ~isempty(str2double(detrndmeth(4))) && ~strcmp(detrndmeth(4),'9')

    % Check a curve-fitting toolbox license is available:
    BF_CheckToolbox('curve_fitting_toolbox');

    [cfun, gof] = fit(r,y,detrndmeth);
    y_fit = feval(cfun,r);
    y_d = y - y_fit;

% 3) Spline detrend
elseif length(detrndmeth) == 8 && strcmp(detrndmeth(1:6),'spline') && ~isempty(str2double(detrndmeth(7))) && ~isempty(str2double(detrndmeth(8)))
    nknots = str2double(detrndmeth(7));
    intp = str2double(detrndmeth(8));

    %% Check that a Curve-Fitting Toolbox license is available:
    BF_CheckToolbox('curve_fitting_toolbox')

	spline = spap2(nknots,intp,r,y); % just a single middle knot with cubic interpolants
	y_spl = fnval(spline,1:N); % evaluate at the 1:N time intervals
    y_d = y - y_spl';

% 4) Differencing
elseif length(detrndmeth) == 5 && strcmp(detrndmeth(1:4),'diff') && ~isempty(str2double(detrndmeth(5)))
    ndiffs = str2double(detrndmeth(5));
    y_d = diff(y,ndiffs); % difference the series n times

% 5) Median Filter
elseif length(detrndmeth) > 7 && strcmp(detrndmeth(1:7),'medianf') && ~isempty(str2double(detrndmeth(8:end)))
    n = str2double(detrndmeth(8:end)); % order of filtering
    y_d = medfilt1(y,n);

% 6) Running Average
elseif length(detrndmeth) > 3 && strcmp(detrndmeth(1:3),'rav') && ~isempty(str2double(detrndmeth(4:end)))
    n = str2double(detrndmeth(4:end)); % the window size
    y_d = filter(ones(1,n)/n,1,y);

% 7) Resample
elseif length(detrndmeth) > 9 && strcmp(detrndmeth(1:9),'resample_')
    % check a valid structure
    ss = textscan(detrndmeth,'%s%n%n','delimiter','_');
    if ~isempty(ss{2}), p = ss{2};
    else return
    end
    if ~isempty(ss{3}), q = ss{3};
    else return
    end
    y_d = resample(y,p,q);

% 8) Log Returns
elseif strcmp(detrndmeth,'logr')
    if all(y > 0), y_d = diff(log(y));
    else
        out = NaN; return % return all NaNs
    end

% 9) Box-Cox Transformation
elseif strcmp(detrndmeth,'boxcox')
    % Requires a financial toolbox to run boxcox, check one is available:
    BF_CheckToolbox('financial_toolbox');

    if all(y > 0), y_d = boxcox(y);
    else
        out = NaN; return % return all NaNs
    end
else
    error('Invalid detrending method ''%s''',detrendmeth)
end

% ------------------------------------------------------------------------------
%% Quick error check
% ------------------------------------------------------------------------------
if all(y_d == 0)
    out = NaN;
    return
end

% ------------------------------------------------------------------------------
%% TESTS ON THE ORIGINAL AND PROCESSED SIGNALS
% ------------------------------------------------------------------------------
% z-score both (these metrics will need it, and not done before-hand
% because of positive-only data, etc.
y = zscore(y);
y_d = zscore(y_d);

% 1) Stationarity

% (a) StatAv
out.statav2 = SY_StatAv(y_d,'seg',2) / SY_StatAv(y,'seg',2);
out.statav4 = SY_StatAv(y_d,'seg',4) / SY_StatAv(y,'seg',4);
out.statav6 = SY_StatAv(y_d,'seg',6) / SY_StatAv(y,'seg',6);
out.statav8 = SY_StatAv(y_d,'seg',8) / SY_StatAv(y,'seg',8);
out.statav10 = SY_StatAv(y_d,'seg',10) / SY_StatAv(y,'seg',10);

% (b) Sliding window mean
out.swms2_2 = SY_SlidingWindow(y_d,'mean','std',2,2) / SY_SlidingWindow(y,'mean','std',2,2);
out.swms5_1 = SY_SlidingWindow(y_d,'mean','std',5,1) / SY_SlidingWindow(y,'mean','std',5,1);
out.swms5_2 = SY_SlidingWindow(y_d,'mean','std',5,2) / SY_SlidingWindow(y,'mean','std',5,2);
out.swms10_1 = SY_SlidingWindow(y_d,'mean','std',10,1) / SY_SlidingWindow(y,'mean','std',10,1);
out.swms10_1 = SY_SlidingWindow(y_d,'mean','std',10,2) / SY_SlidingWindow(y,'mean','std',10,2);

% (c) Sliding window std
out.swss2_1 = SY_SlidingWindow(y_d,'std','std',2,1) / SY_SlidingWindow(y,'std','std',2,1);
out.swss2_2 = SY_SlidingWindow(y_d,'std','std',2,2) / SY_SlidingWindow(y,'std','std',2,2);
out.swss5_1 = SY_SlidingWindow(y_d,'std','std',5,1) / SY_SlidingWindow(y,'std','std',5,1);
out.swss5_2 = SY_SlidingWindow(y_d,'std','std',5,2) / SY_SlidingWindow(y,'std','std',5,2);
out.swss10_1 = SY_SlidingWindow(y_d,'std','std',10,1) / SY_SlidingWindow(y,'std','std',10,1);
out.swss10_2 = SY_SlidingWindow(y_d,'std','std',10,2) / SY_SlidingWindow(y,'std','std',10,2);

% 2) Gaussianity
% (a) kernel density fit
me1 = DN_SimpleFit(y_d,'gauss1',0); % kernel density fit to 1-peak gaussian
me2 = DN_SimpleFit(y,'gauss1',0); % kernel density fit to 1-peak gaussian
if (~isstruct(me1) && isnan(me1)) || (~isstruct(me2) && isnan(me2))
   % fitting gaussian failed -- returns a NaN rather than a structure
    out.gauss1_kd_r2 = NaN;
    out.gauss1_kd_adjr2 = NaN;
    out.gauss1_kd_rmse = NaN;
    out.gauss1_kd_resAC1 = NaN;
    out.gauss1_kd_resAC2 = NaN;
    out.gauss1_kd_resruns = NaN;
else
    out.gauss1_kd_r2 = me1.r2/me2.r2;
    out.gauss1_kd_adjr2 = me1.adjr2/me2.adjr2;
    out.gauss1_kd_rmse = me1.rmse/me2.rmse;
    out.gauss1_kd_resAC1 = me1.resAC1/me2.resAC1;
    out.gauss1_kd_resAC2 = me1.resAC2/me2.resAC2;
    out.gauss1_kd_resruns = me1.resruns/me2.resruns;
end

% (c) compare distribution to fitted normal distribution
me1 = DN_CompareKSFit(y_d,'norm');
me2 = DN_CompareKSFit(y,'norm');

out.kscn_adiff = me1.adiff/me2.adiff;
out.kscn_peaksepy = me1.peaksepy/me2.peaksepy;
out.kscn_peaksepx = me1.peaksepx/me2.peaksepx;
out.kscn_olapint = me1.olapint/me2.olapint;
out.kscn_relent = me1.relent/me2.relent;

% (d) p-values from gof tests
out.htdt_chi2n = HT_DistributionTest(y_d,'chi2gof','norm',10) / HT_DistributionTest(y,'chi2gof','norm',10); % chi2
out.htdt_ksn = HT_DistributionTest(y_d,'ks','norm') / HT_DistributionTest(y,'ks','norm'); % Kolmogorov-Smirnov

% 3) Outliers
out.olbt_m2 = DN_OutlierTest(y_d,2,'mean') / DN_OutlierTest(y,2,'mean');
out.olbt_m5 = DN_OutlierTest(y_d,5,'mean') / DN_OutlierTest(y,5,'mean');
out.olbt_s2 = DN_OutlierTest(y_d,2,'std') / DN_OutlierTest(y,2,'std');
out.olbt_s5 = DN_OutlierTest(y_d,5,'std') / DN_OutlierTest(y,5,'std');

end
