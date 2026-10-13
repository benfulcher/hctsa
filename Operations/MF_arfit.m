function out = MF_arfit(y, pmin, pmax, selector)
% MF_arfit   The coefficients, order, residuals, and oscillation modes of a best-fitting autoregressive model.
%
% Autoregressive (AR) models are fitted with orders p = pmin, pmin + 1, ..., pmax,
% using the ARfit package, with no intercept (the input should be z-scored).
% The optimal model order is selected using Schwarz's Bayesian Criterion (SBC) by
% default; the coefficients, residuals, and modes are those of the model of this
% order. The outputs also give the SBC and the log final prediction error (FPE) at
% every order and how sharp their minima are, summaries of the residuals, 95%
% margins of error on the coefficients, and statistics from an eigendecomposition
% of the fitted model into modes (oscillation periods, damping times, and
% excitations, with times in samples).
%
% ARfit is freely available at http://www.gps.caltech.edu/~tapio/arfit/
%
% ---INPUTS:
% y, the input time series
% pmin, the minimum AR model order to fit (default: 1)
% pmax, the maximum AR model order to fit (default: 10)
% selector, criterion to select the optimal time-series model order ('sbc',
%           the default, or 'fpe'; cf. ARFIT package documentation)
%
% ---OUTPUTS:
% A1, A2, A3, A4, A5, A6: the first six AR coefficients of the selected model (NaN
%       if the selected order is lower)
% maxA, minA, meanA, stdA, sumA: maximum, minimum, mean, standard deviation, and
%       sum of the AR coefficients
% rmsA: square root of the sum of squared AR coefficients (their Euclidean norm)
% C: the estimated noise variance
% sbc_1, sbc_2, ..., sbc_<pmax>: Schwarz's criterion at each order from pmin to pmax
% minsbc, popt_sbc: the minimum SBC and its position within pmin:pmax
% aroundmin_sbc: absolute minimum SBC relative to the mean absolute SBC at the
%       adjacent orders
% fpe_1, fpe_2, fpe_3, fpe_4, fpe_5, fpe_6, fpe_7, fpe_8: log final prediction error at each
%       order (fpe_<k> for order k, k = pmin to pmax)
% minfpe, popt_fpe, aroundmin_fpe: as for the SBC
% res_siglev: p-value of the Li-McLeod portmanteau test of residual autocorrelation
%       (lags up to 20)
% meane, meanabs, stde, maxonstd, ac1, ac2, ac3, propbth, taurat, sws, swm: summaries
%       of the residuals (fit minus data), from MF_ResidualAnalysis at the 'core'
%       level: mean, mean absolute value, standard deviation, largest absolute value
%       in standard deviations, autocorrelation at lags 1 to 3, proportion of the
%       first 25 autocorrelations inside the significance band, ratio of residual to
%       data decorrelation time, and the variability of the residual standard
%       deviation and mean across 5 windows
% aerr_min, aerr_max, aerr_mean: minimum, maximum, and mean 95% margin of error
%       of the AR coefficients
% maxReLambda, maxImLambda, maxabsLambda, stdabsLambda: maximum real part, maximum
%       imaginary part, maximum modulus, and standard deviation of the moduli of the
%       eigenvalues of the companion matrix of the fitted AR model (the roots of its
%       characteristic polynomial; from ARFIT_armode). The largest modulus is the
%       spectral radius, near 1 for a nearly non-stationary process.
% hasInfper: the number of eigenmodes with infinite oscillation period
% meanper, stdper, maxper, minper, meanpererr: mean, standard deviation, maximum,
%       minimum, and mean margin of error of the finite oscillation periods
% meantau, maxtau, mintau, stdtau, meantauerr: mean, maximum, minimum, standard
%       deviation, and mean margin of error of the damping times
% maxexctn, minexctn, meanexctn, stdexctn: maximum, minimum, mean, and standard
%       deviation of the excitations (relative dynamical importance, summing to 1)
%
% ---REFERENCES:
% A. Neumaier and T. Schneider, "Estimation of parameters and eigenmodes of
% multivariate autoregressive models", ACM Trans. Math. Softw. 27, 27 (2001).
% T. Schneider and A. Neumaier, "Algorithm 808: ARFIT---a Matlab package for the
% estimation of parameters and eigenmodes of multivariate autoregressive models",
% ACM Trans. Math. Softw. 27, 58 (2001).
%
% ---NOTES:
% NaN is returned for a (nearly) exactly predictable series: when the estimated noise
% variance is below 1e-12 of the variance of the series (e.g., an exact sinusoid), the
% coefficients are not determined and the residuals are rounding noise.
%
% popt_sbc and popt_fpe are positions within pmin:pmax, so they equal the model
% order only when pmin = 1.

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
%% Check Inputs
% ------------------------------------------------------------------------------
if size(y, 2) > size(y, 1)
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
if ~exist('ARFIT_arfit', 'file')
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

% An exactly predictable series (e.g., a sinusoid) has a singular design and a noise
% variance at the level of rounding error: the coefficients are not determined, the
% residuals are rounding noise, and the eigenmodes sit on the unit circle (infinite
% damping times), so there is nothing to report
if Cest < 1e-12 * var(y) || ~(var(y) > 0)
	out = NaN; return
end

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
% returning the first 6, and NaN if don't exist
out.A1 = Aest(1);
for i = 2:6
	if popt >= i
		out.(sprintf('A%u', i)) = Aest(i);
	else
		% This coefficient was not estimated at the selected order -- report NaN, not 0,
		% which would be indistinguishable from a genuinely-fitted zero coefficient.
		% (Affects ~36% of series at this mop's order range, so this is not an edge case.)
		out.(sprintf('A%u', i)) = NaN;
	end
end

% (ii) Summary statistics on the coefficients
out.maxA = max(Aest);
out.minA = min(Aest);
out.meanA = mean(Aest);
out.stdA = std(Aest);
out.sumA = sum(Aest);
out.rmsA = sqrt(sum(Aest.^2));
% (Dropped: sumsqA = sum(Aest.^2). rmsA is its square root, a strictly monotone transform,
%  so the two are rank-identical by construction.)

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
	out.(sprintf('sbc_%u', ps(i))) = SBC(i);
	% eval(sprintf('out.sbc_%u = SBC(%u);',ps(i),i));
end

% (ii) Return minimum
out.minsbc = min(SBC);
out.popt_sbc = find(SBC == min(SBC), 1, 'first');

% (iii) How convincing is the minimum?
% adjacent values
if (out.popt_sbc > 1) && (out.popt_sbc < length(SBC));
	meanaround = mean(abs([SBC(out.popt_sbc - 1), SBC(out.popt_sbc + 1)]));
elseif out.popt_sbc == 1
	meanaround = abs(SBC(out.popt_sbc + 1)); % just the next value
elseif out.popt_sbc == length(SBC) % really an else
	meanaround = abs(SBC(out.popt_sbc - 1)); % just the previous value
else
	error('Weird error!');
end
out.aroundmin_sbc = abs(min(SBC)) / meanaround;

% ------------------------------------------------------------------------------
% (5) Aikake's Final Prediction Error (FPE)
% ------------------------------------------------------------------------------
% (i) Return all
for i = 1:length(ps)
	out.(['fpe_', num2str(ps(i))]) = FPE(i);
	% eval(sprintf('out.fpe_%u = FPE(%u);',ps(i),i));
end
% (ii) Return minimum
out.minfpe = min(FPE);
out.popt_fpe = find(FPE == min(FPE), 1, 'first');

% (iii) How convincing is the minimum?
% adjacent values
if out.popt_fpe > 1 && out.popt_fpe < length(FPE);
	meanaround = mean(abs([FPE(out.popt_fpe - 1), FPE(out.popt_fpe + 1)]));
elseif out.popt_fpe == 1
	meanaround = abs(FPE(out.popt_fpe + 1)); % just the next value
elseif out.popt_fpe == length(FPE) % really an else
	meanaround = abs(FPE(out.popt_fpe - 1));
else
	error('Weird error!!');
end
out.aroundmin_fpe = abs(min(FPE)) / meanaround;

% -------------------------------------------------------------------------------
%% (II) Test Residuals
% -------------------------------------------------------------------------------

% Run code from ARfit package:
[siglev, res] = ARFIT_arres(west, Aest, y);

% (1) Significance Level
out.res_siglev = siglev;

% (2) Residual diagnostics, through the shared contract at the cheap 'core' level.
% Replaces the hand-rolled res_ac1, res_ac1_norm and pcorr_res: the first two are now
% ac1 (ac1n was dropped as ~92% recoverable from the rest), and pcorr_res -- the
% proportion of the first 20 autocorrelations exceeding 1.96/sqrt(N) -- is superseded by
% propbth, which is the same idea over 25 lags at the 2.6/sqrt(N) threshold.
residOut = MF_ResidualAnalysis(-res, y, 'core'); % ARFIT_arres returns data minus fit; the contract is prediction minus data
fields = fieldnames(residOut);
for k = 1:length(fields)
	out.(fields{k}) = residOut.(fields{k});
end


% -------------------------------------------------------------------------------
%% (III) Confidence Intervals
% -------------------------------------------------------------------------------

% Run code from ARfit package:
Aerr = ARFIT_arconf(Aest, Cest, th);

% Return mean/min/max error margins
out.aerr_min = min(Aerr);
out.aerr_max = max(Aerr);
out.aerr_mean = mean(Aerr);

% -------------------------------------------------------------------------------
%% (III) Eigendecomposition
% -------------------------------------------------------------------------------

% Run code from the ARfit package
[~, ~, per, tau, exctn, lambda] = ARFIT_armode(Aest, Cest, th);

% lambda: eigenvalues of the companion matrix of the AR model (complex in conjugate
%         pairs for oscillatory modes; modulus < 1 for a stable model)
% Serr: +/- margins of error (95% confidence intervals)
% per: periods of oscillation (margins of error in second row)
% tau: damping times (margins of error in second row)
% exct: measures of relative dynamical importance of eigenmodes

% Since there will be a variable number, best to just use summaries
% (These were previously computed from S, the last component of each unit-length,
%  phase-adjusted eigenvector, which reflects the eigenvector normalization rather than
%  the dynamics; they are now computed from the eigenvalues, lambda, themselves.)
out.maxReLambda = max(real(lambda));
out.maxImLambda = max(imag(lambda));
out.maxabsLambda = max(abs(lambda));
out.stdabsLambda = std(abs(lambda));

% Often you get infinite periods of oscillation -- remove these for the purposes
% of taking stats:
perSpecial = ~isfinite(per);
perFiltered = per;
perFiltered(perSpecial) = NaN;

out.hasInfper = sum(perSpecial(1, :));
out.meanper = mean(perFiltered(1, :),'omitnan');
out.stdper = std(perFiltered(1, :),0,'omitnan');
% (max and min ignore NaNs by default)
out.maxper = max(perFiltered(1, :));
out.minper = min(perFiltered(1, :));
out.meanpererr = mean(per(2, :),'omitnan');

out.meantau = mean(tau(1, :));
out.maxtau = max(tau(1, :));
out.mintau = min(tau(1, :));
out.stdtau = std(tau(1, :));
out.meantauerr = mean(tau(2, :));

out.maxexctn = max(exctn);
out.minexctn = min(exctn);
out.meanexctn = mean(exctn);
out.stdexctn = std(exctn);

end
