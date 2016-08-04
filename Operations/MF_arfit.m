function out = MF_arfit(y,pmin,pmax,selector)
% MF_arfit      Statistics of a fitted AR model to a time series.
%
% Uses various functions implemented in the ARfit package, which is
% freely-available at http://www.gps.caltech.edu/~tapio/arfit/
%
% cf. "Estimation of parameters and eigenmodes of multivariate autoregressive
%       models", A. Neumaier and T. Schneider, ACM Trans. Math. Softw. 27, 27 (2001)
%
% cf. "Algorithm 808: ARFIT---a Matlab package for the estimation of parameters
%      and eigenmodes of multivariate autoregressive models",
%      T. Schneider and A. Neumaier, ACM Trans. Math. Softw. 27, 58 (2001)
%
% Autoregressive (AR) models are fitted with orders p = pmin, pmin + 1, ..., pmax.
%
% The optimal model order is selected using Schwartz's Bayesian Criterion (SBC).
%
%---INPUTS:
% y, the input time series
% pmin, the minimum AR model order to fit
% pmax, the maximum AR model order to fit
% selector, crierion to select optimal time-series model order (e.g., 'sbc', cf.
%           ARFIT package documentation)
%
%---OUTPUTS: include the model coefficients obtained, the SBCs at each model
% order, various tests on residuals, and statistics from an eigendecomposition
% of the time series using the estimated AR model.

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
%% Check Inputs
% ------------------------------------------------------------------------------
if size(y,2) > size(y,1)
    y = y'; % needs to be a column vector
end
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
if ~exist('ARFIT_arfit','file')
    error('Cannot find the function ''ARFIT_arfit''. There''s a problem with the ARfit toolbox.')
end

% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------
%% (I) Fit AR model
% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------
% Run the code with no intercept vector (all input data should be
% zero-mean, z-scored)
[west, Aest, Cest, SBC, FPE, th] = ARFIT_arfit(y, pmin, pmax, selector, 'zero');

% First, some definitions
ps = (pmin:pmax);
popt = length(Aest);

% ------------------------------------------------------------------------------
% (1) Intercept west
% ------------------------------------------------------------------------------
% west = 0 -- as specified

% ------------------------------------------------------------------------------
% (2) Coefficients Aest
% ------------------------------------------------------------------------------
% (i) Return the raw coefficients
% somewhat problematic since will depend on order fitted. We can try
% returning the first 6, and 0s if don't exist
out.A1 = Aest(1);
for i = 2:6
    if popt >= i
        out.(sprintf('A%u',i)) = Aest(i);
    else
        % it's as if the higher order coefficients are all zero
        out.(sprintf('A%u',i)) = 0;
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

% ------------------------------------------------------------------------------
% (3) Noise covariance matrix, Cest
% ------------------------------------------------------------------------------
% In our case of a univariate time series, just a scalar for the noise
% magnitude.
out.C = Cest;

% ------------------------------------------------------------------------------
% (4) Schwartz's Bayesian Criterion, SBC
% ------------------------------------------------------------------------------
% (not included in default HCTSA library -- rather the FPE is used)
% There will be a value for each model order from pmin:pmax
% (i) Return all
for i = 1:length(ps)
    out.(sprintf('sbc_%u',ps(i))) = SBC(i);
    % eval(sprintf('out.sbc_%u = SBC(%u);',ps(i),i));
end

% (ii) Return minimum
out.minsbc = min(SBC);
out.popt_sbc = find(SBC == min(SBC),1,'first');

% (iii) How convincing is the minimum?
% adjacent values
if (out.popt_sbc > 1) && (out.popt_sbc < length(SBC));
    meanaround = mean(abs([SBC(out.popt_sbc-1), SBC(out.popt_sbc+1)]));
elseif out.popt_sbc == 1
    meanaround = abs(SBC(out.popt_sbc+1)); % just the next value
elseif out.popt_sbc == length(SBC) % really an else
    meanaround = abs(SBC(out.popt_sbc-1)); % just the previous value
else
    error('Weird error!');
end
out.aroundmin_sbc = abs(min(SBC))/meanaround;

% ------------------------------------------------------------------------------
% (5) Aikake's Final Prediction Error (FPE)
% ------------------------------------------------------------------------------
% (i) Return all
for i = 1:length(ps)
    out.(['fpe_',num2str(ps(i))]) = FPE(i);
    % eval(sprintf('out.fpe_%u = FPE(%u);',ps(i),i));
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

%-------------------------------------------------------------------------------
%% (II) Test Residuals
%-------------------------------------------------------------------------------

% Run code from ARfit package:
[siglev, res] = ARFIT_arres(west,Aest,y);

% (1) Significance Level
out.res_siglev = siglev;

% (2) Correlation test of residuals
% error margins are within 1.96/sqrt(N);
out.res_ac1 = CO_AutoCorr(res,1,'Fourier');
out.res_ac1_norm = CO_AutoCorr(res,1,'Fourier')/sqrt(N); % normalize by sqrt(N)

% Calculate correlations up to 20, return how many exceed significance threshold
acf = CO_AutoCorr(res,1:20,'Fourier');
out.pcorr_res = sum(abs(acf)>1.96/sqrt(N))/20;

%-------------------------------------------------------------------------------
%% (III) Confidence Intervals
%-------------------------------------------------------------------------------

% Run code from ARfit package:
Aerr = ARFIT_arconf(Aest, Cest, th);

% Return mean/min/max error margins
out.aerr_min = min(Aerr);
out.aerr_max = max(Aerr);
out.aerr_mean = mean(Aerr);

%-------------------------------------------------------------------------------
%% (III) Eigendecomposition
%-------------------------------------------------------------------------------

% Run code from the ARfit package
[S, ~, per, tau, exctn] = ARFIT_armode(Aest, Cest, th);

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

% Often you get infinite periods of oscillation -- remove these for the purposes
% of taking stats:
perSpecial = ~isfinite(per);
perFiltered = per;
perFiltered(perSpecial) = NaN;

out.hasInfper = sum(perSpecial(1,:));
out.meanper = nanmean(perFiltered(1,:));
out.stdper = nanstd(perFiltered(1,:));
out.maxper = nanmax(perFiltered(1,:));
out.minper = nanmin(perFiltered(1,:));
out.meanpererr = nanmean(per(2,:));

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
