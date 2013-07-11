function out = MF_arfit(y,pmin,pmax,selector)
% Uses MATLAB code from ARfit to fit an AR model and return statistics about it.
% http://www.gps.caltech.edu/~tapio/arfit/
% Adapted by Ben Fulcher 28/1/2010

%% Check Inputs
N = length(y); % time series length

if nargin < 2 || isempty(pmin)
    pmin = 1;
end
if nargin < 3 || isempty(pmax)
    pmax = 10;
end
if nargin < 4 || isempty(selector)
     selector = 'sbc';
     % Use Schwartz's Bayesian Criterion to choose optimum model order
end

% Check the ARfit toolbox is installed and in the Matlab path
myarfit = which('arfit');
if isempty(myarfit)
    error('Cannot find the function ''arfit''. Have you installed the ARFIT toolbox?')
end

%% (I) Fit AR model
% Run the code with no intercept vector (all input data should be
% zero-mean, z-scored)
[west, Aest, Cest, SBC, FPE, th] = arfit(y, pmin, pmax, selector, 'zero');


% (0) First, some definitions
ps = (pmin:pmax);
popt = length(Aest);

% (1) Intercept west
% west = 0 -- as specified


% (2) Coefficients Aest
% (i) Return the raw coefficients
% somewhat problematic since will depend on order fitted. We can try
% returning the first 6, and 0s if don't exist
out.A1 = Aest(1);
for i = 2:6
    if popt >= i
        eval(sprintf('out.A%u = Aest(%u);',i,i))
    else
        % it's as if the higher order coefficients are all zero
        eval(sprintf('out.A%u = 0;',i,i))
    end
end

% (ii) Summary statistics on the coefficients
out.maxA = max(Aest);
out.minA = min(Aest);
out.meanA = mean(Aest);
out.stdA = std(Aest);
out.sumA = sum(Aest);
out.rmsA = sqrt(sum(Aest.^2));
out.sumsqA = sum(Aest.^2);

% (3) Noise covariance matrix, Cest
% In our case of a univariate time series, just a scalar for the noise
% magnitude.
out.C = Cest;

% (4) Schwartz's Bayesian Criterion, SBC
% There will be a value for each model order from pmin:pmax
% (i) Return all
for i = 1:length(ps)
    eval(sprintf('out.sbc_%u = SBC(%u);',ps(i),i));
end

% (ii) Return minimum
out.minsbc = min(SBC);
out.popt_sbc = find(SBC == min(SBC),1,'first');

% (iii) How convincing is the minimum?
% adjacent values
if out.popt_sbc > 1 && out.popt_sbc < length(SBC);
    meanaround = mean(abs([SBC(out.popt_sbc-1),SBC(out.popt_sbc+1)]));
elseif out.popt_sbc == 1
    meanaround = abs(SBC(out.popt_sbc+1)); % just the next value
elseif out.popt_sbc == length(SBC) % really an else
    meanaround = abs(SBC(out.popt_sbc-1)); % just the previous value
else
    error('Weird error!');
end
out.aroundmin_sbc = abs(min(SBC))/meanaround;


% (5) Aikake's Final Prediction Error (FPE)
% (i) Return all
for i = 1:length(ps)
    eval(sprintf('out.fpe_%u = FPE(%u);',ps(i),i));
end
% (ii) Return minimum
out.minfpe = min(FPE);
out.popt_fpe = find(FPE == min(FPE),1,'first');

% (iii) How convincing is the minimum?
% adjacent values
if out.popt_fpe > 1 && out.popt_fpe < length(FPE);
    meanaround = mean(abs([FPE(out.popt_fpe-1),FPE(out.popt_fpe+1)]));
elseif out.popt_fpe == 1
    meanaround = abs(FPE(out.popt_fpe+1)); % just the next value
elseif out.popt_fpe == length(FPE) % really an else
    meanaround = abs(FPE(out.popt_fpe-1));
else
    error('Weird error!!');
end
out.aroundmin_fpe = abs(min(FPE))/meanaround;


%% (II) Test Residuals

% Run code from ARfit package:
[siglev, res] = arres(west,Aest,y);

% (1) Significance Level
out.res_siglev = siglev;

% (2) Correlation test of residuals
% error margins are within 1.96/sqrt(N);
out.res_ac1 = CO_autocorr(res,1);
out.res_ac1_norm = CO_autocorr(res,1)/sqrt(N); % normalize by sqrt(N)

% Calculate correlations up to 20, return how many exceed significance threshold
acf = CO_autocorr(res,1:20);
out.pcorr_res = sum(abs(acf)>1.96/sqrt(N))/20;

%% (III) Confidence Intervals

% Run code from ARfit package:
Aerr = arconf(Aest, Cest, th);

% Return mean/min/max error margins
out.aerr_min = min(Aerr);
out.aerr_max = max(Aerr);
out.aerr_mean = mean(Aerr);

%% (III) Eigendecomposition

% Run code from the ARfit package
[S, Serr, per, tau, exctn] = armode(Aest, Cest, th);

% S: eigenmodes
% Serr: +/- margins of error (95% confidence intervals)
% per: periods of oscillation (margins of error in second row)
% tau: damping times (margins of error in second row)
% exct: measures of relative dynamical importance of eigenmodes

% Since there will be a variable number, best to just use summaries
out.maxReS = max(real(S));
out.maxImS = max(imag(S));
out.maxabsS = max(abs(S));
out.stdabsS = std(abs(S));

out.meanper = mean(per(1,:));
out.stdper = std(per(1,:));
out.maxper = max(per(1,:));
out.minper = min(per(1,:));
out.meanpererr = mean(per(2,:));

out.meantau = mean(tau(1,:));
out.maxtau = max(tau(1,:));
out.mintau = min(tau(1,:));
out.stdtau = std(tau(1,:));
out.meantauerr = mean(tau(2,:));

out.maxexctn = max(exctn);
out.minexctn = min(exctn);
out.meanexctn = mean(exctn);
out.stdexctn = std(exctn);

end