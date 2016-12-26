function out = SP_Summaries(y,psdmeth,wmeth,nf,dologabs)
% SP_Summaries  Statistics of the power spectrum of a time series
%
% The estimation can be done using a periodogram, using the periodogram code in
% Matlab's Signal Processing Toolbox, or a fast fourier transform, implemented
% using Matlab's fft code.
%
%---INPUTS:
% y, the input time series
%
% psdmeth, the method of obtaining the spectrum from the signal:
%               (i) 'periodogram': periodogram
%               (ii) 'fft': fast fourier transform
%
% wmeth, the window to use:
%               (i) 'boxcar'
%               (iii) 'bartlett'
%               (iv) 'hann'
%               (v) 'hamming'
%               (vi) 'none'
%
% nf, the number of frequency components to include, if
%           empty (default), it's approx length(y)
%
% dologabs, if 1, takes log amplitude of the signal before
%           transforming to the frequency domain.
%
% doPower, analyzes the power spectrum rather than amplitudes of a Fourier
%          transform
%
%---OUTPUTS:
% Statistics summarizing various properties of the spectrum,
% including its maximum, minimum, spread, correlation, centroid, area in certain
% (normalized) frequency bands, moments of the spectrum, Shannon spectral
% entropy, a spectral flatness measure, power-law fits, and the number of
% crossings of the spectrum at various amplitude thresholds.

% ------------------------------------------------------------------------------
% Copyright (C) 2016, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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
%% Check that a Curve-Fitting Toolbox license is available:
% ------------------------------------------------------------------------------
BF_CheckToolbox('curve_fitting_toolbox')

% ------------------------------------------------------------------------------
% Check inputs, set defaults:
% ------------------------------------------------------------------------------
if size(y,2) > size(y,1);
    y = y'; % Time series must be a column vector
end
if nargin < 2 || isempty(psdmeth)
    psdmeth = 'fft'; % fft by default
end
if nargin < 3 || isempty(wmeth)
    wmeth = 'hamming'; % Hamming window by default
end
if nargin < 4
    nf = [];
end
if nargin < 5 || isempty(dologabs)
    dologabs = 0;
end

if dologabs % a boolean
    % Analyze the spectrum of logarithmic absolute deviations
    y = log(abs(y));
end

doPlot = 0; % plot outputs
Ny = length(y); % time-series length

%-------------------------------------------------------------------------------
% Set window (for periodogram and welch):
%-------------------------------------------------------------------------------
if ismember(psdmeth,{'periodogram','welch'})
    switch wmeth % method to use for the window
        case 'none'
            window = [];
        case 'hamming'
            window = hamming(Ny);
        case 'hann'
            window = hann(Ny);
        case 'bartlett'
            window = bartlett(Ny);
        case 'boxcar'
            window = boxcar(Ny);
        case 'rect'
            window = rectwin(Ny);
        otherwise
            % There are other options, but these aren't implemented here
            error('Unknown window ''%s''',wmeth);
    end
end

% ------------------------------------------------------------------------------
% Compute the Fourier Transform
% ------------------------------------------------------------------------------
switch psdmeth
    case 'periodogram'
        if isempty(nf)
            % (2) Estimate the spectrum
            [S, w] = periodogram(y,window);
        else
            w = linspace(0,pi,nf);
            [S, w] = periodogram(y,window,w);
        end

    case 'fft'
        % Fast Fourier Transform
        Fs = 1; % sampling frequency
        NFFT = 2^nextpow2(Ny);
        f = Fs/2*linspace(0,1,NFFT/2+1); % frequency
        w = 2*pi*f'; % angular frequency (as column vector)
        S = fft(y,NFFT); % Fourier Transform
        S = 2*abs(S(1:NFFT/2+1)).^2/Ny; % single-sided power spectral density
        S = S/(2*pi); % convert to angular frequency space

    case 'welch'
        % Welch power spectral density estimate:
        Fs = 1; % sampling frequency
        N = 2^nextpow2(Ny);
        [S, f] = pwelch(y,window,[],N,Fs);
        w = 2*pi*f'; % angular frequency
        S = S/(2*pi); % adjust so that area remains normalized in angular frequency space

    otherwise
        error('Unknown spectral estimation method ''%s''',psdmeth);
end

if ~any(isfinite(S)) % no finite values in the power spectrum
    % This time series must be really weird -- return NaN (unsuitable operation)...
    fprintf(1,'NaN in power spectrum? A weird time series.\n');
    out = NaN; return
end

% Ensure both w and S are row vectors:
if size(S,1) > size(S,2)
    S = S';
end
if size(w,1) > size(w,2)
    w = w';
end

if doPlot
    figure('color','w')
    plot(w,S,'.-k'); % plot the spectrum
    % Area under S should sum to 1 if a power spectral density estimate:
    title(sprintf('Area under psd curve = %.1f (= %.1f)',sum(S*(w(2)-w(1))),var(y)));
end

N = length(S); % = length(w)
logS = log(S);
dw = w(2) - w(1); % spacing increment in w

% ------------------------------------------------------------------------------
% Simple measures of the power spectrum
% ------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Peaks:
%-------------------------------------------------------------------------------
% Maximum, and max peak width:
[out.maxS, i_maxS] = max(S);
out.maxw = w(i_maxS);
out.maxWidth = w(i_maxS + find(S(i_maxS+1:end) < out.maxS,1,'first')) - ...
                    w(find(S(1:i_maxS-1) < out.maxS,1,'last'));
if isempty(out.maxWidth);
    out.maxWidth = 0;
end

% Characterize all peaks using findpeaks function:
% Minimum angular separation of 0.02...?
minDist_w = 0.02;
ptsPerw = length(S)/pi;
minPkDist = ceil(minDist_w*ptsPerw);
[pkHeight,pkLoc,pkWidth,pkProm] = findpeaks(S,'SortStr','descend','minPeakDistance',minPkDist);
pkWidth = pkWidth/ptsPerw;
pkLoc = pkLoc/ptsPerw;

% Characterize mean peak prominence
% (use prominence threshold of 2...?)
out.numPeaks = length(pkHeight); % total number of peaks
out.numPromPeaks_1 = sum(pkProm > 1); % number of peaks with prominence of at least 1
out.numPromPeaks_2 = sum(pkProm > 2); % number of peaks with prominence of at least 2
out.numPromPeaks_5 = sum(pkProm > 5); % number of peaks with prominence of at least 5
out.numPeaks_overmean = sum(pkProm>mean(pkProm)); % number of peaks with prominence greater than the mean (low for skewed distn)
out.maxProm = max(pkProm); % maximum prominence of any peak
out.meanProm_2 = mean(pkProm(pkProm > 2)); % mean peak prominence of those with prominence of at least 2

out.meanPeakWidth_prom2 = mean(pkWidth(pkProm > 2)); % mean peak width of peaks with prominence of at least 2
out.width_weighted_prom = sum(pkWidth.*pkProm)/sum(pkProm);

% Power in top N peaks:
nn = @(x) 1:min(x,out.numPeaks);
out.peakPower_2 = sum(pkHeight(nn(2)).*pkWidth(nn(2)));
out.peakPower_5 = sum(pkHeight(nn(5)).*pkWidth(nn(5)));
out.peakPower_prom2 = sum(pkHeight(pkProm > 2).*pkWidth(pkProm > 2)); % power in peaks with prominence of at least 2
out.w_weighted_peak_prom = sum(pkLoc.*pkProm)/sum(pkProm); % where are prominent peaks located on average (weighted by prominence)
out.w_weighted_peak_height = sum(pkLoc.*pkHeight)/sum(pkHeight); % where are prominent peaks located on average (weighted by height)

% Number of peaks required to get to 50% of power in peaks
peakPower = pkHeight.*pkWidth;
out.numPeaks_50power = find(cumsum(peakPower) > 0.5*sum(peakPower),1,'first');
out.peakpower_1 = peakPower(1)/sum(peakPower);

%-------------------------------------------------------------------------------
% Distribution
%-------------------------------------------------------------------------------
% Quantiles:
out.iqr = iqr(S);
out.logiqr = iqr(logS);
out.q25 = quantile(S,0.25);
out.median = median(S);
out.q75 = quantile(S,0.75);

% Moments:
out.std = std(S);
out.stdlog = log(out.std);
out.logstd = std(logS);
out.mean = mean(S);
out.logmean = mean(logS);
for i = 3:5
    out.(sprintf('mom%u',i)) = DN_Moments(S,i);
end

% Autocorr:
autoCorrs_S = CO_AutoCorr(S,1:4,'Fourier');
out.ac1 = autoCorrs_S(1);
out.ac2 = autoCorrs_S(2);
out.tau = CO_FirstZero(S,'ac');

% ------------------------------------------------------------------------------
% Shape of cumulative sum curve
% ------------------------------------------------------------------------------
csS = cumsum(S);

f_frac_w_max = @(f) w(find(csS >= csS(end)*f,1,'first'));

% At what frequency is csS a fraction p of its maximum?
out.wmax_5 = f_frac_w_max(0.05);
out.wmax_10 = f_frac_w_max(0.1);
out.wmax_25 = f_frac_w_max(0.25);
out.centroid = f_frac_w_max(0.5);
out.wmax_75 = f_frac_w_max(0.75);
out.wmax_90 = f_frac_w_max(0.9);
out.wmax_95 = f_frac_w_max(0.95);
out.wmax_99 = f_frac_w_max(0.99);

% Width of saturation measures
out.w10_90 = out.wmax_90 - out.wmax_10; % from 10% to 90%:
out.w25_75 = out.wmax_75 - out.wmax_25;

% ------------------------------------------------------------------------------
% Fit some functions to this cumulative sum:
% ------------------------------------------------------------------------------
% (i) Quadratic
[c, gof] = fit(w',csS','poly2');
out.fpoly2csS_p1 = c.p1;
out.fpoly2csS_p2 = c.p2;
out.fpoly2csS_p3 = c.p3;
out.fpoly2_sse = gof.sse;
out.fpoly2_r2 = gof.rsquare;
out.fpoly2_rmse = gof.rmse;

% (ii) Fit polysat a*x^2/(b+x^2) (has zero derivative at zero, though)
s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[csS(end), 100]);
f = fittype('a*x^2/(b+x^2)','independent','x','options',s); % set 'a' from maximum
[c, gof] = fit(w',csS',f);
out.fpolysat_a = c.a;
out.fpolysat_b = c.b; % this is important
out.fpolysat_r2 = gof.rsquare; % this is more important!
out.fpolysat_rmse = gof.rmse;

% ------------------------------------------------------------------------------
% Shannon spectral entropy
% ------------------------------------------------------------------------------
Hshann = -S.*log(S); % Shannon function
out.spect_shann_ent = sum(Hshann);
out.spect_shann_ent_norm = mean(Hshann);

% ------------------------------------------------------------------------------
% "Spectral Flatness Measure"
% ------------------------------------------------------------------------------
% which is given in dB as 10 log_10(gm/am) where gm is the geometric mean and am
% is the arithmetic mean of the power spectral density
out.sfm = 10*log10(geomean(S)/mean(S));

% ------------------------------------------------------------------------------
% Areas under power spectrum
% ------------------------------------------------------------------------------
% Area up to peak: (may be more appropriate in squared log units?)
out.areatopeak = sum(S(1:i_maxS))*dw;
out.ylogareatopeak = sum(logS(1:i_maxS))*dw; % (semilogy)
% out.logareatopeak=sum(logS(1:i_maxS).*diff(logw(1:i_maxS+1)));

% ------------------------------------------------------------------------------
%% Robust (e.g., iteratively re-weighted least squares) linear fits to log-log
%   plots
% ------------------------------------------------------------------------------

% Use the local function, giveMeRobustStats, which adds robust
% stats to the out structure.

% Suppress rank deficient warnings for this section:
warning('off','stats:robustfit:RankDeficient')

% (1): Across full range
r_all = (w > 0); % avoid -Inf for log(0) when w = 0;
out = giveMeRobustStats(log(w(r_all)),log(S(r_all)),'linfitloglog_all',out);

% (2): First half (low frequency)
r_lf = (w > 0); % w(1) = 0 -> log(0) = -Inf
r_lf(floor(N/2)+1:end) = 0; % remove second half of angular frequenciesf
out = giveMeRobustStats(log(w(r_lf)),log(S(r_lf)),'linfitloglog_lf',out);

% (3): Second half (high frequency)
r_hf = floor(N/2)+1:N;
out = giveMeRobustStats(log(w(r_hf)),log(S(r_hf)),'linfitloglog_hf',out);

% (4): Middle half (mid-frequencies)
r_mf = round(N/4):round(N*3/4);
out = giveMeRobustStats(log(w(r_mf)),log(S(r_mf)),'linfitloglog_mf',out);

% (5) Fit linear to semilog plot (across full range)
out = giveMeRobustStats(w,log(S),'linfitsemilog_all',out);

% Turn the rank-deficient warnings back on
warning('on','stats:robustfit:RankDeficient')

% ------------------------------------------------------------------------------
%% Power in specific frequency bands
% ------------------------------------------------------------------------------
% *** DO THIS BY BUFFER COMMAND: AND WHILE AT IT LOOK AT STATIONARITY OF
% POWER SPECTRUM....

% 2 bands
split = buffer(S,floor(N/2));
if size(split,2) > 2, split = split(:,1:2); end
out.area_2_1 = sum(split(:,1))*dw;
out.logarea_2_1 = sum(log(split(:,1)))*dw;
out.area_2_2 = sum(split(:,2))*dw;
out.logarea_2_2 = sum(log(split(:,2)))*dw;
out.statav2_m = std(mean(split))/std(S);
out.statav2_s = std(std(split))/std(S);

% 3 bands
split = buffer(S,floor(N/3));
if size(split,2) > 3, split = split(:,1:3); end
out.area_3_1 = sum(split(:,1))*dw;
out.logarea_3_1 = sum(log(split(:,1)))*dw;
out.area_3_2 = sum(split(:,2))*dw;
out.logarea_3_2 = sum(log(split(:,2)))*dw;
out.area_3_3 = sum(split(:,3))*dw;
out.logarea_3_3 = sum(log(split(:,3)))*dw;
out.statav3_m = std(mean(split))/std(S);
out.statav3_s = std(std(split))/std(S);

% 4 bands
split = buffer(S,floor(N/4));
if size(split,2) > 4, split = split(:,1:4); end
out.area_4_1 = sum(split(:,1))*dw;
out.logarea_4_1 = sum(log(split(:,1)))*dw;
out.area_4_2 = sum(split(:,2))*dw;
out.logarea_4_2 = sum(log(split(:,2)))*dw;
out.area_4_3 = sum(split(:,3))*dw;
out.logarea_4_3 = sum(log(split(:,3)))*dw;
out.area_4_4 = sum(split(:,4))*dw;
out.logarea_4_4 = sum(log(split(:,4)))*dw;
out.statav4_m = std(mean(split))/std(S);
out.statav4_s = std(std(split))/std(S);

% 5 bands
split = buffer(S,floor(N/5));
if size(split,2) > 5, split = split(:,1:5); end
out.area_5_1 = sum(split(:,1))*dw;
out.logarea_5_1 = sum(log(split(:,1)))*dw;
out.area_5_2 = sum(split(:,2))*dw;
out.logarea_5_2 = sum(log(split(:,2)))*dw;
out.area_5_3 = sum(split(:,3))*dw;
out.logarea_5_3 = sum(log(split(:,3)))*dw;
out.area_5_4 = sum(split(:,4))*dw;
out.logarea_5_4 = sum(log(split(:,4)))*dw;
out.area_5_5 = sum(split(:,5))*dw;
out.logarea_5_5 = sum(log(split(:,5)))*dw;
out.statav5_m = std(mean(split))/std(S);
out.statav5_s = std(std(split))/std(S);

% ------------------------------------------------------------------------------
% Count crossings:
% Get a horizontal line and count the number of crossings with the power spectrum
% ------------------------------------------------------------------------------
ncrossfn_rel = @(f) sum(BF_sgnchange(S - f*max(S)));

out.ncross_f05 = ncrossfn_rel(0.05);
out.ncross_f01 = ncrossfn_rel(0.1);
out.ncross_f02 = ncrossfn_rel(0.2);
out.ncross_f05 = ncrossfn_rel(0.5);

%-------------------------------------------------------------------------------
% function mel = w2mel(w) % convert to mel spectrum
%     mel = 1127*log(w/(1400*pi)+1);
% end

function out = giveMeRobustStats(xData,yData,textID,out)
    % Add statistics to the output structure from a robust linear fit
    % between xData and yData

    % Perform the fit:
    [a, stats] = robustfit(xData,yData);

    % Add the statistics to the output structure:
    out.(sprintf('%s_a1',textID)) = a(1); % robust intercept
    out.(sprintf('%s_a2',textID)) = a(2); % robust gradient
    % ratio of sigma estimates between ordinary least squares (ols) and the robust fit:
    out.(sprintf('%s_sigrat',textID)) = stats.ols_s/stats.robust_s;
    % esimate of sigma as the larger of robust_s and a weighted average of ols_s and robust_s:
    out.(sprintf('%s_sigma',textID)) = stats.s;
    out.(sprintf('%s_sea1',textID)) = stats.se(1); % standard error of 1st coefficient estimate
    out.(sprintf('%s_sea2',textID)) = stats.se(2); % standard error of 2nd coefficient estimate
end

end
