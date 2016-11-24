function out = SD_TSTL_surrogates(y, tau, nsurr, surrMethod, surrfn, randomSeed)
% SD_TSTL_surrogates    Surrogate time-series analysis
%
% Generates surrogate time series and tests them against the original time
% series according to some test statistics: T_{C3}, using the TSTOOL code tc3 or
% T_{rev}, using TSTOOL code trev.
%
%---INPUTS:
% y, the input time series
%
% tau, the autocorrelation lag length <x_n x_{n-tau} x_{n-2tau)>/abs(<x_n
%                                                   x_{n-tau}|^3/2
% nsurr, the number of surrogates to generate
%
% surrmeth, the method of generating surrogates:
%               (i) 1: randomizes phases of fourier spectrum
%               (ii) 2:  (see Theiler algorithm II)
%               (iii) 3: permutes samples randomly
%
% surrfn, the surrogate statistic to evaluate on all surrogates, either 'tc3' or
%           'trev'
%
% randomSeed, whether (and how) to reset the random seed, using BF_ResetSeed
%
%---OUTPUTS: include the Gaussianity of the test statistics, a z-test, and
% various tests based on fitted kernel densities.
%
% TSTOOL: http://www.physik3.gwdg.de/tstool/

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
%% Check inputs, set defaults
% ------------------------------------------------------------------------------
% 1) time delay, TAU
if nargin < 2 || isempty(tau)
    tau = 1;
end
if strcmp(tau,'ac')
    tau = CO_FirstZero(y,'ac');
elseif strcmp(tau,'mi')
    tau = CO_FirstMin(y,'mi');
end

% 2) number of surrogate data sets to generate, NSURR
if nargin < 3 || isempty(nsurr)
    nsurr = 50;
end

% 3) surrogate data method, SURRMETHOD
if nargin < 4 || isempty(surrMethod)
    fprintf(1,'Surrogate method set to default: ''surrogate1''.\n');
    surrMethod = 1;
end
% surrMethod = 1: randomizes phases of fourier spectrum
% surrMethod = 2:  (see Theiler algorithm II)
% surrMethod = 3: permutes samples randomly

% 4) surrogate function, SURRFN
if nargin < 5 || isempty(surrfn)
    fprintf(1,'surrogate function set to default value: ''tc3''.\n');
    surrfn = 'tc3';
end

% 5) randomSeed: how to treat the randomization
if nargin < 6
    randomSeed = []; % default for BF_ResetSeed
end

% ------------------------------------------------------------------------------
%% Do the calculation
% ------------------------------------------------------------------------------
% Make a TSTOOL signal object of time series
s = signal(y);

% Control random seed (for reproducibility):
BF_ResetSeed(randomSeed);

switch surrfn
    case 'tc3'
        % Run external TSTOOL code, tc3
        rs = tc3(s, tau, nsurr, surrMethod);
    case 'trev'
        % Run external TSTOOL code, trev
        rs = trev(s, tau, nsurr, surrMethod);
    otherwise
        error('Unknown surrogate function ''%s''',surrfn)
end

tc3dat = data(rs);
if all(isnan(tc3dat))
    error('TSTOOL: ''%s'' failed',surrfn);
end
tc3_y = tc3dat(1);
tc3_surr = tc3dat(2:end);

% figure; histogram(tc3_surr);

% ------------------------------------------------------------------------------
%% Get some outputs
% ------------------------------------------------------------------------------
% These are completely from my intuition

% 1) fit a Gaussian to surrogates
% [muhat,sigmahat] = normfit(tc3_surr);
muhat = mean(tc3_surr);
sigmahat = std(tc3_surr);
% probability of data given Guassian surrogates
out.normpatponmax = normpdf(tc3_y,muhat,sigmahat)/normpdf(muhat,muhat,sigmahat);

% Probability at least that distance from mean, using ztest:

% 2) stds from mean
out.stdfrommean = abs(tc3_y - mean(tc3_surr))/std(tc3_surr);
% (~equivalent to a z-test:)
[~, out.ztestp] = ztest(tc3_y, muhat, sigmahat);
% (both of these stats are a monotonic function of normpatponmax)

% iqrs from median
out.iqrsfrommedian = abs(tc3_y - median(tc3_surr))/iqr(tc3_surr);

% 3) basic info on surrogates
out.stdsurr = sigmahat;
out.meansurr = muhat;

% 4) kernel density test
[ksf, ksx] = ksdensity(tc3_surr,'function','pdf');
% hold on;plot(ksx,ksf,'r')
ksdx = ksx(2) - ksx(1);
ihit = find(ksx > tc3_y,1,'first');

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
