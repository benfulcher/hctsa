% SP_Summaries
% 
% Returns a set of measures summarizing an estimate of the Fourier transform of
% the signal.
% 
% The estimation can be done using a periodogram, using the periodogram code in
% Matlab's Signal Processing Toolbox, or a fast fourier transform, implemented
% using Matlab's fft code.
% 
% INPUTS:
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
% dopower, analyzes the power spectrum rather than amplitudes of a Fourier
%          transform
% 
% Outputs are statistics summarizing various properties of the spectrum,
% including its maximum, minimum, spread, correlation, centroid, area in certain
% (normalized) frequency bands, moments of the spectrum, Shannon spectral
% entropy, a spectral flatness measure, power-law fits, and the number of
% crossings of the spectrum at various amplitude thresholds.
% 
% ------------------------------------------------------------------------------
% Copyright (C) 2013,  Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones., "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2010). DOI: 10.1098/rsif.2013.0048
%
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
% details.
% 
% You should have received a copy of the GNU General Public License along with
% this program.  If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

function out = SP_Summaries(y,psdmeth,wmeth,nf,dologabs,dopower)
% Ben Fulcher, August 2009

%% Check that a Curve-Fitting Toolbox license is available:
BF_CheckToolbox('curve_fitting_toolbox')

% Check inputs, set defaults:
if size(y,2) > size(y,1); y = y'; end % time series must be a column vector
if nargin < 2 || isempty(psdmeth), psdmeth = 'fft'; end
if nargin < 3 || isempty(wmeth), wmeth = 'none'; end % no window by default
if nargin < 4, nf = []; end
if nargin < 5 || isempty(dologabs), dologabs = 0; end
if nargin < 6 || isempty(dopower), dopower = 1; end

if dologabs % a boolean
    y = log(abs(y));
    % This analyzes is the spectrum of logarithmic absolute deviations
end

Ny = length(y); % time-series length

switch psdmeth
    case 'periodogram'
        % (1) set the window
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
            otherwise
                % There are other options, but these aren't implemented here
                error('Unknown window ''%s''',wmeth);
        end
        
        if isempty(nf)
            % (2) Estimate the spectrum
            [S, w] = periodogram(y,window);
        else
            w = linspace(0,pi,nf);
            [S, w] = periodogram(y,window,w);
        end
        
    case 'fft'
        Fs = 1; % sampling frequency
        NFFT = 2^nextpow2(Ny);
        f = Fs/2*linspace(0,1,NFFT/2+1); % frequency
        w = 2*pi*f'; % angular frequency as column vector
        S = fft(y,NFFT)/Ny; % Fourier Transform
        S = 2*abs(S(1:NFFT/2+1)); % single-sided amplitudes
        % convert to power spectrum later if necessary

    otherwise
        error('Unknwon power spectral density method ''%s''',psdmeth);
end

if ~any(isfinite(S)) % no finite values in the power spectrum
    % This time series must be really weird -- return NaN (unsuitable operation)...
    fprintf(1,'This is a weird time series\n');
    out = NaN; return
end

% Look at power spectrum rather than amplitudes
if dopower
    S = S.^2;
end

if size(S,1) > size(S,2); S = S'; w = w'; end
N = length(S); % = length(w)
% plot(w,S); % plot the spectrum
logS = log(S);
logw = log(w);
dw = w(2) - w(1);
% Normalize to 1: necessary if input not z-scored
% S=S/(sum(S)*dw);

% Simple measures
out.std = std(S);
out.stdlog = log(out.std);
out.logstd = std(logS);
out.mean = mean(S);
out.meanlog = log(out.mean);
out.logmean = mean(logS);
[out.maxS, i_maxS] = max(S);
out.maxSlog = log(out.maxS);
out.maxw = w(i_maxS);
out.median = median(S);
out.medianlog = log(out.median);
out.melmax = w2mel(out.median); % transform to mel scale
out.iqr = iqr(S);
out.logiqr = iqr(logS);
out.q25 = quantile(S,0.25);
out.q25log = log(out.q25);
out.q75 = quantile(S,0.75);
out.q75log = log(out.q75);
out.ac1 = CO_AutoCorr(S,1);
out.ac2 = CO_AutoCorr(S,1);
out.ac3 = CO_AutoCorr(S,1);
out.ac4 = CO_AutoCorr(S,1);
out.tau = CO_FirstZero(S,'ac');
out.logac1 = CO_AutoCorr(logS,1);
out.logac2 = CO_AutoCorr(logS,1);
out.logac3 = CO_AutoCorr(logS,1);
out.logac4 = CO_AutoCorr(logS,1);
out.logtau = CO_FirstZero(logS,'ac');
out.logmaxonlogmean1e = out.maxSlog/mean(logS(logS<out.maxSlog/exp(1)));
out.maxwidth = w(i_maxS+find(logS(i_maxS+1:end)<out.maxSlog/exp(1),1,'first'))-...
                w(find(logS(1:i_maxS)<out.maxSlog/exp(1),1,'last'));
if isempty(out.maxwidth); out.maxwidth=0; end
% out.maxlogwidth=log(w(find(logS(i_maxS+1:end)<out.maxSlog/exp(1),1,'first')))-...
%                 log(w(find(logS(1:i_maxS-1)<out.maxSlog/exp(1),1,'last')));

% input('see cum sum?')
csS = cumsum(S);
% plot(w,csS);

% Measures of central location
out.centroid = w(find(csS > csS(end)/2,1,'first')); % where area under curve is same above
                                            % and below this frequency

% Shape of cumulative sum curve
% 1) Quantiles
% where is csS at fraction p of its maximum?
out.q1 = w(find(csS > 0.01*csS(end),1,'first'));
out.q1mel = w2mel(out.q1);
out.q5 = w(find(csS > 0.05*csS(end),1,'first'));
out.q5mel = w2mel(out.q5);
out.q10 = w(find(csS > 0.10*csS(end),1,'first'));
out.q10mel = w2mel(out.q10);
out.q25 = w(find(csS > 0.25*csS(end),1,'first'));
out.q25mel = w2mel(out.q25);
out.q50 = w(find(csS > 0.50*csS(end),1,'first')); % centroid
out.q50mel = w2mel(out.q50);
out.q75 = w(find(csS > 0.75*csS(end),1,'first'));
out.q75mel = w2mel(out.q75);
out.q90 = w(find(csS > 0.90*csS(end),1,'first'));
out.q90mel = w2mel(out.q90);
out.q95 = w(find(csS > 0.95*csS(end),1,'first'));
out.q95mel = w2mel(out.q95);
out.q99 = w(find(csS > 0.99*csS(end),1,'first'));
out.q99mel = w2mel(out.q99);

% width of saturation measures
out.w1_99 = out.q99-out.q1;
out.w1_99mel = w2mel(out.w1_99);
out.w5_95 = out.q95-out.q5;
out.w5_95mel = w2mel(out.w5_95);
out.w10_90 = out.q90-out.q10;% from 10% to 90%:
out.w10_90mel = w2mel(out.w10_90);
out.w25_75 = out.q75-out.q25;
out.w25_75mel = w2mel(out.w25_75);

% 2) Moments of the power spectrum
% take moments from 3--9
for i = 3:9
    themom = DN_Moments(S,i);
    eval(sprintf('out.mom%u = themom;',i));
end

% fit some functions to this cumulative sum:
% (i) quadratic
[c, gof] = fit(w',csS','poly2');
out.fpoly2csS_p1 = c.p1;
out.fpoly2csS_p2 = c.p2;
out.fpoly2csS_p3 = c.p3;
out.fpoly2_sse = gof.sse;
out.fpoly2_r2 = gof.rsquare;
out.fpoly2_adjr2 = gof.adjrsquare;
out.fpoly2_rmse = gof.rmse;

% (ii) fit polysat ~x^2/(b+x^2) (has zero derivative at zero, though)
% s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[1]);
% f = fittype('a*x^2/(b+x^2)','problem','a','independent','x','options',s); % set 'a' from maximum
% [c,gof] = fit(w,csS,f,'problem',csS(end)); % the saturation value
s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[csS(end), 100]);
f = fittype('a*x^2/(b+x^2)','independent','x','options',s); % set 'a' from maximum
[c, gof] = fit(w',csS',f);
out.fpolysat_a = c.a;
out.fpolysat_b = c.b; % this is important
out.fpolysat_r2 = gof.rsquare; % this is more important!
out.fpolysat_adjr2 = gof.adjrsquare;
out.fpolysat_rmse = gof.rmse;



% KP on http://www.dsprelated.com/showmessage/108326/1.php
% If zscored, power spectrum should already be normalized
% sum(P(f))(dw)=1  (where P(f) is the power spectrum)
% 
% 2) Transform with the Shannon function:
% H(f)=Q(f)[log(1/Q(f))]
% 
% 3) Spectral entropy:
% E=sum(H(f))/log(Nf)  (where Nf is the number of frequency components.
% 
% For wavelet spectral entropy, your step (1) is appropriate.  You need to 
% replace you step (2) with:
% 
% 2) Transform with the Shannon function: 
% Hi=Pi[log(1/Pi)]
% 
% 3) Wavelet spectral entropy:
% E=sum(Hi)/log(Ni) (where Ni is the number of subbands)

% Shannon spectral entropy
Hshann = S.*log(1./S); % Shannon function
out.spect_shann_ent = sum(Hshann);
out.spect_shann_ent_norm = sum(Hshann)/length(S);

% "Spectral Flatness Measure" which is given in dB
% as 10 log_10(gm/am) where gm is the geometric mean and am is the arithmetic
% mean of the power spectrum.
out.sfm = 10*log10(geomean(S)/mean(S));


% Areas under power spectrum
% Area up to peak: (may be more appropriate in squared log units?)
out.areatopeak = sum(S(1:i_maxS))*dw;
out.ylogareatopeak = sum(logS(1:i_maxS))*dw; % (semilogy)
% out.logareatopeak=sum(logS(1:i_maxS).*diff(logw(1:i_maxS+1)));


%% robust linear fits to log-log plots
% (1): full range

warning('off','stats:robustfit:RankDeficient') % suppress these warnings

[a, stats] = robustfit(log(w),log(S));
out.linfitloglog_all_a1 = a(1); % robust intercept
out.linfitloglog_all_a2 = a(2); % robust gradient
out.linfitloglog_all_sigrat = stats.ols_s/stats.robust_s;
out.linfitloglog_all_s = stats.s; % esimate on sigma
out.linfitloglog_all_sea1 = stats.se(1); % standard error of coefficient 1 estimate
out.linfitloglog_all_sea2 = stats.se(2); % standard error of coefficient 2 estimate

% (2): first half (low frequency)
[a, stats] = robustfit(log(w(1:round(N/2))),log(S(1:round(N/2))));
out.linfitloglog_lf_a1 = a(1); % robust intercept
out.linfitloglog_lf_a2 = a(2); % robust gradient
out.linfitloglog_lf_sigrat = stats.ols_s/stats.robust_s;
out.linfitloglog_lf_s = stats.s;
out.linfitloglog_lf_sea1 = stats.se(1);
out.linfitloglog_lf_sea2 = stats.se(2);

% (3): second half (high frequency)
[a, stats] = robustfit(log(w(round(N/2):end)),log(S(round(N/2):end)));
out.linfitloglog_hf_a1 = a(1); % robust intercept
out.linfitloglog_hf_a2 = a(2); % robust gradient
out.linfitloglog_hf_sigrat = stats.ols_s/stats.robust_s;
out.linfitloglog_hf_s = stats.s;
out.linfitloglog_hf_sea1 = stats.se(1);
out.linfitloglog_hf_sea2 = stats.se(2);

% (4): middle half (mid-frequencies)
[a, stats] = robustfit(log(w(round(N/4):round(N*3/4))),log(S(round(N/4):round(N*3/4))));
out.linfitloglog_mf_a1 = a(1); % robust intercept
out.linfitloglog_mf_a2 = a(2); % robust gradient
out.linfitloglog_mf_sigrat = stats.ols_s/stats.robust_s;
out.linfitloglog_mf_s = stats.s;
out.linfitloglog_mf_sea1 = stats.se(1);
out.linfitloglog_mf_sea2 = stats.se(2);

% fit linear to semilog plot (full range)
[a, stats] = robustfit(w,log(S));
out.linfitloglog_all_a1 = a(1); % robust intercpet
out.linfitloglog_all_a2 = a(2); % robust gradient
out.linfitloglog_all_sigrat = stats.ols_s/stats.robust_s;
out.linfitloglog_all_s = stats.s;
out.linfitloglog_all_sea1 = stats.se(1);
out.linfitloglog_all_sea2 = stats.se(2);

warning('on','stats:robustfit:RankDeficient') % turn these warnings back on


%% Power in bands
% *** DO THIS BY BUFFER COMMAND: AND WHILE AT IT LOOK AT STATIONARITY OF
% POWER SPECTRUM....
% 2 bands
split = buffer(S,floor(N/2));
if size(split,2) > 2, split = split(:,1:2); end
out.area_2_1 = sum(split(:,1))*dw;
out.logarea_2_1 = sum(log(split(:,1)))*dw;
out.area_2_2 = sum(split(:,2))*dw;
out.logarea_2_2 = sum(log(split(:,2)))*dw;
out.rawstatav2_m = std(mean(split));
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
out.rawstatav3_m = std(mean(split));
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
out.rawstatav4_m = std(mean(split));
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
out.rawstatav5_m = std(mean(split));
out.statav5_m = std(mean(split))/std(S);
out.statav5_s = std(std(split))/std(S);

% How many crossings
% Give a horizontal line and count the # crossings with the power spectrum
ncrossfn = @(x) sum(BF_sgnchange(S - x));

out.ncross01 = ncrossfn(0.1);
out.ncross02 = ncrossfn(0.2);
out.ncross05 = ncrossfn(0.5);
out.ncross1 = ncrossfn(1);
out.ncross2 = ncrossfn(2);
out.ncross5 = ncrossfn(5);
out.ncross10 = ncrossfn(10);
out.ncross20 = ncrossfn(20);

function mel = w2mel(w) % convert to mel
    mel = 1127*log(w/(1400*pi)+1);
end

end