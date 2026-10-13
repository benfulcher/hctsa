function out = MF_armax(y, orders, pTrain, numSteps)
% MF_armax   The coefficients of a fitted ARMA model, and how well it predicts the later part of the series.
%
% Fits an autoregressive moving-average (ARMA) model with orders [p, q] to the whole
% time series by conditional least squares from a Hannan-Rissanen start (BF_ARMAfit). The
% coefficients, their uncertainties, and the goodness of fit are from this fit.
% The model is then fitted again to the first pTrain proportion of the time series
% and used to predict the remainder numSteps samples ahead; the prediction
% residuals (prediction minus data) are summarized with MF_ResidualAnalysis.
%
% ---INPUTS:
%
% y, the input time series
%
% orders, a two-vector for p and q, the AR and MA components of the model,
%           respectively (default: [3, 3])
%
% pTrain, the proportion of data to train the model on (the remainder is used
%           for testing; default: 0.8)
%
% numSteps, number of steps to predict into the future for testing the model
%           (default: 1)
%
% ---OUTPUTS:
% From the model fitted to the entire time series, in the convention
% y(t) + a1 y(t-1) + ... + ap y(t-p) = e(t) + c1 e(t-1) + ... + cq e(t-q):
% AR_1, AR_2, AR_3: the AR coefficients a1, ..., ap (the negatives of the usual AR
%       coefficients)
% MA_1, MA_2: the MA coefficients c1, ..., cq
% maxda, maxdc: the largest standard error of the AR and MA coefficients (from the
%       Jacobian of the residuals at the optimum; NaN if p or q is 0)
% noisevar, lossfn, fpe: the noise variance lossfn*N/(N - p - q), the loss function
%       (the mean squared one-step prediction error), and Akaike's final prediction
%       error lossfn*(1 + (p + q)/N)/(1 - (p + q)/N) of the fit
% From the residuals of the predictions of the held-out portion (MF_ResidualAnalysis):
% meane, meanabs, stde, maxonstd: mean, mean absolute value, standard deviation,
%       and largest absolute value (in standard deviations) of the residuals
% ac1, ac2, ac3: residual autocorrelation at lags 1 to 3
% propbth: proportion of the first 25 residual autocorrelations within the
%       significance band (|r| < 2.6/sqrt(N))
% ftbth: the first lag at which the residual autocorrelation is within that band
% taurat: ratio of the residual decorrelation time to the data decorrelation time
% sws, swm: variability of the residual standard deviation and mean across 5 windows
% normksstat: Kolmogorov-Smirnov statistic of the residuals against a Gaussian
% popt, minsbc: the order of the best AR model fitted to the residuals (chosen by
%       SBC, from 1 to 10) and its SBC
%
% ---NOTES:
% The fit (BF_ARMAfit) is a fully specified deterministic procedure: a Hannan-Rissanen
% start refined by Levenberg-Marquardt iterations on the conditional sum of squares,
% with the MA polynomial kept invertible. The coefficients are poorly determined where an ARMA model is poorly
% identified (nearly cancelling AR and MA polynomials, as for near-white series, where
% only the loss, noisevar and fpe are well determined).
% Predictions start from a diffuse (completely unknown) initial state (BF_SSPredict), so
% the first max(p, q) + numSteps - 1 samples of a series are not predicted; the loss and
% the residual summaries average over the samples predicted.
% The held-out portion starts at sample floor(pTrain*N), overlapping the training
% portion by one sample. All outputs are NaN if the model cannot be fitted (a series
% too short for the orders).

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
%% Check inputs
% ------------------------------------------------------------------------------
y = y(:);
N = length(y); % number of samples

if nargin < 2 || isempty(orders)
	orders = [3, 3]; % AR3, MA3
end
p = orders(1); q = orders(2);
if nargin < 3 || isempty(pTrain)
	pTrain = 0.8; % train on 80% of the data
end
if nargin < 4 || isempty(numSteps)
	numSteps = 1; % one-step-ahead predictions
end

% ------------------------------------------------------------------------------
%% Fit the ARMA model to the whole time series
% ------------------------------------------------------------------------------
try
	[a, c, coefSD, loss] = BF_ARMAfit(y, p, q);
catch
	out = NaN; return % too few samples for the orders
end

for i = 1:p
	out.(sprintf('AR_%u', i)) = a(i + 1);
end
for i = 1:q
	out.(sprintf('MA_%u', i)) = c(i + 1);
end
if p == 0
	out.maxda = NaN;
else
	out.maxda = max(coefSD(1:p));
end
if q == 0
	out.maxdc = NaN;
else
	out.maxdc = max(coefSD(p+1:end));
end

% Goodness of fit
nPar = p + q;
out.noisevar = loss * N / (N - nPar);
out.lossfn = loss;
out.fpe = loss * (1 + nPar/N) / (1 - nPar/N); % Akaike's final prediction error

% ------------------------------------------------------------------------------
%% Fit to the training portion, predict the test portion
% ------------------------------------------------------------------------------
ytrain = y(1:floor(pTrain * N));
ytest = y(floor(pTrain * N):end); % overlap by one sample
try
	[ap, cp] = BF_ARMAfit(ytrain, p, q);
catch
	out = NaN; return
end
[A, K, C] = BF_ARMAtoSS(ap, cp);
yp = BF_SSPredict(A, K, C, ytest, numSteps);
mresiduals = yp - ytest; % prediction minus data (the MF_ResidualAnalysis convention)
mresiduals = mresiduals(~isnan(mresiduals)); % the samples that can be predicted

% ------------------------------------------------------------------------------
%% Residual analysis
% ------------------------------------------------------------------------------
residout = MF_ResidualAnalysis(mresiduals, ytest, 'full');
fields = fieldnames(residout);
for k = 1:length(fields)
	out.(fields{k}) = residout.(fields{k});
end

end
