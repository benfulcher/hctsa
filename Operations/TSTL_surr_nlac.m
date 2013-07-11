function out = TSTL_surr_nlac(y, tau, nsurr, surrmethod, surrfn)
% Uses TSTOOL code tc3/trev to generate surrogates and comparing to nonlinear
% autocorrelation statistic
% INPUTS:
% tau: the autocorrelation lag length <x_n x_{n-tau} x_{n-2tau)>/abs(<x_n
% x_{n-tau}|^3/2
% Ben Fulcher 15/11/2009

%% Check inputs, set defaults
% 1) time delay, TAU
if nargin < 2 || isempty(tau)
    tau = 1;
end
if strcmp(tau,'ac')
    tau = CO_fzcac(y);
elseif strcmp(tau,'mi')
    tau = CO_fmmi(y);
end


% 2) number of surrogate data sets to generate, NSURR
if nargin < 3 || isempty(nsurr)
    nsurr = 50;
end

% 3) surrogate data method, SURRMETHOD
if nargin < 4 || isempty(surrmethod)
    disp('you should set the surrogate method.');
    disp('just this once I''ll do it for you -- surrogate1');
    surrmethod = 1;
end
% surrmethod = 1: randomizes phases of fourier spectrum
% surrmethod = 2:  (see Theiler algorithm II)
% surrmethod = 3: permutes samples randomly

% 4) surrogate function, SURRFN
if nargin < 5 || isempty(surrfn)
    disp('Using tc3 by default');
    surrfn = 'tc3';
end

%% Do the calculation
% Make a signal object of time series
s = signal(y);

switch surrfn
    case 'tc3'
        rs = tc3(s, tau, nsurr, surrmethod);
    case 'trev'
        rs = trev(s, tau, nsurr, surrmethod);
end

tc3dat = data(rs);
if all(isnan(tc3dat))
    disp('Failed horribly');
    out = NaN; return
end
tc3_y = tc3dat(1);
tc3_surr = tc3dat(2:end);

[n, x] = hist(tc3_surr,50);
% hold off; plot(x,n); hold on; plot(tc3_y,max(n),'or'); hold off;
    
%% Get some outputs
% these are completely made up by me

% 1) fit a Gaussian to surrogates
% [muhat,sigmahat] = normfit(tc3_surr);
muhat = mean(tc3_surr);
sigmahat = std(tc3_surr);
% probability of data given Guassian surrogates
% out.normpatp = normpdf(tc3_y,muhat,sigmahat);
out.normpatponmax = normpdf(tc3_y,muhat,sigmahat)/normpdf(muhat,muhat,sigmahat);
% probability at least that distance from mean:
[ztest_h ztest_p] = ztest(tc3_y, muhat, sigmahat);
out.ztestp = ztest_p;

% 2) stds from mean
out.stdfrommean = abs(tc3_y - mean(tc3_surr))/std(tc3_surr);
% iqrs from median
out.iqrsfrommedian = abs(tc3_y - median(tc3_surr))/iqr(tc3_surr);

% 3) basic info on surrogates
out.stdsurr = sigmahat;
out.meansurr = muhat;

% 4) kernel density test
% ksx = linspace(min(tc3_surr),max(tc3_surr),200); % don't know how to pick
% extremeties this way...
% ksf = ksdensity(tc3_surr,ksx,'function','pdf');
[ksf ksx] = ksdensity(tc3_surr,'function','pdf');
% hold on;plot(ksx,ksf,'r')
ksdx = ksx(2)-ksx(1);
ihit = find(ksx>tc3_y,1,'first');

if isempty(ihit) %% off the scale!
    out.kspminfromext = 0;
    out.ksphereonmax = 0;
else % on the scale!
    pfromleft = ksdx*sum(ksf(1:ihit));
    % pfromright = ksdx*sum(ksf(ihit+1:end))
    out.kspminfromext = min([pfromleft 1-pfromleft]);
    % out.phereonstd = ksf(ihit)/sigmahat;
    out.ksphereonmax = ksf(ihit)/normpdf(muhat,muhat,sigmahat);
%     out.ksiqrsfrommode = abs(ksx(imode)-ksx(ihit))/iqr(tc3_surr);
end

% iqrs from mode
imode = find(ksf == max(ksf),1,'first');
out.ksiqrsfrommode = abs(ksx(imode)-tc3_y)/iqr(tc3_surr);

end