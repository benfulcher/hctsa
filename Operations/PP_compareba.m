function out=PP_compareba(y,detrndmeth)
% Inputs: y: the time series
%         detrndmeth: the detrending method(s)
% Coding begun on 9/7/09 by Ben Fulcher

%% List of valid detrndmeth:
% If multiple methods, should be in a cell of strings; the methods will be
% executed in order, e.g., {'poly1','sin1'} does a linear polynomial then a
% simple one-frequency seasonal detrend (in that order)

% Later we should take some smarter approaches that choose the
% preprocessing based on some measures. Here we just state what to do and
% do it.

% 1) Polynomials
%   (a) fit polynomial of given order
%       'poly1', 'poly2', 'poly3', 'poly4', 'poly5', 'poly6', 'poly7', 'poly8', 'poly9'
%   (b) fit best polynomial
%       'polybest': determines 'best' by various tests (e.g., whiteness of
%               residuals, etc.)
%   (c) 'fitstrong' only fits if a 'strong' trend

% 2) Sines
%   (a) fit a sine series of a given order
% Fits a form like: a1*sin(b1*x+c1)+a2*sin(b2*x+c2)+...
% Additional number determines how many terms to include in the series:
%   'sin1', 'sin2', 'sin3', 'sin4', 'sin5', 'sin6', 'sin7', 'sin8'
%   (b) fit only if a strong trend (i.e., if the amplitudes are above a
%                                           given threshold)
%   'sin_st1'

% 3) Splines
% Of form 'spline<nknots><interpolant_order>'
% e.g., 'spline45' uses four knots and 5th order interpolants
% (Implements a least squares spline via the spline toolbox function spap2)

% 4) Differencing
% (a)   Take just a simple differencing of the time series
%       of form 'diff<ndiff>' where ndiff is the number of differencings to
%       perform. e.g., 'diff3' performs three recursive differences

% 5) Median filter
% Just take a running median
% Of form 'medianf<n>' where n is the window length.
% e.g., 'medianf3' takes a running median using the median of every 3
%       consecutive values as a point in the filtered time series
% Uses the signal processing toolbox function medfilt1

% 6) Running mean
% Uses the filter function in matlab to perform a running average of order
% n. Of form 'rav<n>' where n is the order of the running average.

% 7) Resample
% Resamples the data using the resample function in MATLAB.
% Of form 'resample_<p>_<q>', where the ratio p/q is the new sampling rate
% e.g., 'resample_1_2' will downsample the signal by one half
% e.g., resample_10_1' will resample the signal to 10 times its original
%                       length

% 8) Log returns
% 'logr'
% Only valid for positive data; otherwise returns a NaN

% 9) Box-Cox transformation
% 'boxcox'
% Only valid for positive only data; otherwise returns a NaN

%% FOREPLAY
N = length(y); % length of time series
% y=zscore(y); % No -- apply transformation to *raw* data
r = (1:N)'; % the x-range over which to fit

%% PREPROCESSINGS:
% DETRENDINGS:
% Do the detrending; converting from y (raw) to y_d (detrended) by
% subtracting some fit y_fit

% 1) Polynomial detrend
% starts with 'poly' and ends with integer from 1--9
if length(detrndmeth) == 5 && strcmp(detrndmeth(1:4),'poly') && ~isempty(str2double(detrndmeth(5)))
    [cfun,gof] = fit(r,y,detrndmeth);
    y_fit = feval(cfun,r);
    y_d = y-y_fit;


% 2) Seasonal detrend
elseif length(detrndmeth) == 4 && strcmp(detrndmeth(1:3),'sin') && ~isempty(str2double(detrndmeth(4))) && ~strcmp(detrndmeth(4),'9')
    [cfun,gof] = fit(r,y,detrndmeth);
    y_fit = feval(cfun,r);
    y_d = y-y_fit;

% 3) Spline detrend
elseif length(detrndmeth) == 8 && strcmp(detrndmeth(1:6),'spline') && ~isempty(str2double(detrndmeth(7))) && ~isempty(str2double(detrndmeth(8)))
    nknots = str2double(detrndmeth(7));
    intp = str2double(detrndmeth(8));
    try
		spline = spap2(nknots,intp,r,y); % just a single middle knot with cubic interpolants
    catch me
		disp(me.message)
		return
	end
	y_spl = fnval(spline,1:N); % evaluate at the 1:N time intervals
    y_d = y-y_spl';
    
% 4) Differencing
elseif length(detrndmeth) == 5 && strcmp(detrndmeth(1:4),'diff') && ~isempty(str2double(detrndmeth(5)))
    ndiffs = str2double(detrndmeth(5));
    y_d = diff(y,ndiffs); % difference the series n times

% 5) Median Filter
elseif length(detrndmeth)>7 && strcmp(detrndmeth(1:7),'medianf') && ~isempty(str2double(detrndmeth(8:end)))
    n = str2double(detrndmeth(8:end)); % order of filtering
    y_d = medfilt1(y,n);

% 6) Running Average
elseif length(detrndmeth)>3 && strcmp(detrndmeth(1:3),'rav') && ~isempty(str2double(detrndmeth(4:end)))
    n = str2double(detrndmeth(4:end)); % the window size
    y_d = filter(ones(1,n)/n,1,y);

% 7) Resample
elseif length(detrndmeth)>9 && strcmp(detrndmeth(1:9),'resample_')
    % check a valid structure
    ss = textscan(detrndmeth,'%s%n%n','delimiter','_');
    if ~isempty(ss{2}), p=ss{2};
    else return
    end
    if ~isempty(ss{3}), q=ss{3};
    else return
    end
    y_d = resample(y,p,q);

% 8) Log Returns
elseif strcmp(detrndmeth,'logr')
    if all(y>0), y_d=diff(log(y));
    else
        out = NaN; % return all NaNs
        return
    end

% 9) Box-Cox Transformation
elseif strcmp(detrndmeth,'boxcox')
    if all(y>0), y_d=boxcox(y);
    else
        out = NaN; % return all NaNs
        return
    end
else
    disp('Invalid detrending method'); return
end

%% quick error checks
if all(y_d == 0)
    out = NaN; return
end

%% TESTS ON THE ORIGINAL AND PROCESSED SIGNALS
% z-score both (these metrics will need it, and not done before-hand
% because of positive-only data, etc.
y = zscore(y); y_d = zscore(y_d);

% 1) Stationarity

% (a) StatAv
out.statav2 = SY_StatAv(y_d,2,'seg')/SY_StatAv(y,2,'seg');
out.statav4 = SY_StatAv(y_d,4,'seg')/SY_StatAv(y,4,'seg');
out.statav6 = SY_StatAv(y_d,6,'seg')/SY_StatAv(y,6,'seg');
out.statav8 = SY_StatAv(y_d,8,'seg')/SY_StatAv(y,8,'seg');
out.statav10 = SY_StatAv(y_d,10,'seg')/SY_StatAv(y,10,'seg');

% (b) Sliding window mean
out.swms2_1 = SY_slidwin(y_d,'mean','std',2,1)/SY_slidwin(y,'mean','std',2,1);
out.swms2_2 = SY_slidwin(y_d,'mean','std',2,2)/SY_slidwin(y,'mean','std',2,2);
out.swms5_1 = SY_slidwin(y_d,'mean','std',5,1)/SY_slidwin(y,'mean','std',5,1);
out.swms5_2 = SY_slidwin(y_d,'mean','std',5,2)/SY_slidwin(y,'mean','std',5,2);
out.swms10_1 = SY_slidwin(y_d,'mean','std',10,1)/SY_slidwin(y,'mean','std',10,1);
out.swms10_1 = SY_slidwin(y_d,'mean','std',10,2)/SY_slidwin(y,'mean','std',10,2);

% (c) Sliding window std
out.swss2_1 = SY_slidwin(y_d,'std','std',2,1)/SY_slidwin(y,'std','std',2,1);
out.swss2_2 = SY_slidwin(y_d,'std','std',2,2)/SY_slidwin(y,'std','std',2,2);
out.swss5_1 = SY_slidwin(y_d,'std','std',5,1)/SY_slidwin(y,'std','std',5,1);
out.swss5_2 = SY_slidwin(y_d,'std','std',5,2)/SY_slidwin(y,'std','std',5,2);
out.swss10_1 = SY_slidwin(y_d,'std','std',10,1)/SY_slidwin(y,'std','std',10,1);
out.swss10_2 = SY_slidwin(y_d,'std','std',10,2)/SY_slidwin(y,'std','std',10,2);

% 2) Gaussianity
% (a) kernel density fit
me1 = MF_M_mtlbfit(y_d,'gauss1',0); % kernel density fit to 1-peak gaussian
me2 = MF_M_mtlbfit(y,'gauss1',0); % kernel density fit to 1-peak gaussian
if (~isstruct(me1) && isnan(me1)) || (~isstruct(me2) && isnan(me2))
   % fitting gaussian failed -- returns a NaN rather than a structure
    out.gauss1_kd_r2 = NaN;
    out.gauss1_kd_adjr2 = NaN;
    out.gauss1_kd_rmse = NaN;
    out.gauss1_kd_resAC1 = NaN;
    out.gauss1_kd_resAC2 = NaN;
    out.gauss1_kd_resruns = NaN;
    % out.gauss1_kd_reslbq=me1.reslbq/me2.reslbq;
else
    out.gauss1_kd_r2 = me1.r2/me2.r2;
    out.gauss1_kd_adjr2 = me1.adjr2/me2.adjr2;
    out.gauss1_kd_rmse = me1.rmse/me2.rmse;
    out.gauss1_kd_resAC1 = me1.resAC1/me2.resAC1;
    out.gauss1_kd_resAC2 = me1.resAC2/me2.resAC2;
    out.gauss1_kd_resruns = me1.resruns/me2.resruns;
    % out.gauss1_kd_reslbq=me1.reslbq/me2.reslbq;
end

%   (b) histogram 10 bins
me1 = MF_M_mtlbfit(y_d,'gauss1',10); % 10-bin histogram fit to 1-peak gaussian
me2 = MF_M_mtlbfit(y,'gauss1',10); % 10-bin histogram fit to 1-peak gaussian
if (~isstruct(me1) && isnan(me1)) || (~isstruct(me2) && isnan(me2))
    out.gauss1_h10_r2 = NaN;
    out.gauss1_h10_adjr2 = NaN;
    out.gauss1_h10_rmse = NaN;
    out.gauss1_h10_resAC1 = NaN;
    out.gauss1_h10_resAC2 = NaN;
    out.gauss1_h10_resruns = NaN;
%     out.gauss1_h10_reslbq = NaN;
else
    out.gauss1_h10_r2 = me1.r2/me2.r2;
    out.gauss1_h10_adjr2 = me1.adjr2/me2.adjr2;
    out.gauss1_h10_rmse = me1.rmse/me2.rmse;
    out.gauss1_h10_resAC1 = me1.resAC1/me2.resAC1;
    out.gauss1_h10_resAC2 = me2.resAC2/me2.resAC2;
    out.gauss1_h10_resruns = me1.resruns/me2.resruns;
%     out.gauss1_h10_reslbq = me1.reslbq/me2.reslbq;
end

% (c) compare distribution to fitted normal distribution
me1 = DN_M_kscomp(y_d,'norm');
me2 = DN_M_kscomp(y,'norm');

out.kscn_adiff = me1.adiff/me2.adiff;
out.kscn_peaksepy = me1.peaksepy/me2.peaksepy;
out.kscn_peaksepx = me1.peaksepx/me2.peaksepx;
out.kscn_olapint = me1.olapint/me2.olapint;
out.kscn_relent = me1.relent/me2.relent;

% (d) p-values from gof tests
out.htdt_chi2n = HT_disttests(y_d,'chi2gof','norm',10)/HT_disttests(y,'chi2gof','norm',10); % chi2
out.htdt_ksn = HT_disttests(y_d,'ks','norm')/HT_disttests(y,'ks','norm'); % Kolmogorov-Smirnov
out.htdt_llfn = HT_disttests(y_d,'lillie','norm')/HT_disttests(y,'ks','norm');

% 3) Outliers
out.olbt_m2 = OL_bentest(y_d,2,1)/OL_bentest(y,2,1);
out.olbt_m5 = OL_bentest(y_d,5,1)/OL_bentest(y,5,1);
out.olbt_s2 = OL_bentest(y_d,2,2)/OL_bentest(y,2,2);
out.olbt_s5 = OL_bentest(y_d,5,2)/OL_bentest(y,5,2);



end
