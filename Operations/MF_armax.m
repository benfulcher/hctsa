function out = MF_armax(y, orders, pTrain, numSteps)
% MF_armax   The coefficients of a fitted ARMA model, and how well it predicts the later part of the series.
%
% Fits an autoregressive moving-average (ARMA) model with orders [p, q] to the whole
% time series, using armax from Matlab's System Identification Toolbox. The
% coefficients, their uncertainties, and the goodness of fit are from this fit.
% The model is then fitted again to the first pTrain proportion of the time series
% and used to predict the remainder numSteps samples ahead; the prediction
% residuals (prediction minus data) are summarized with MF_ResidualAnalysis.
%
% Uses the functions iddata, armax, aic, and predict from Matlab's System
% Identification Toolbox
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
% From the model fitted to the entire time series, in the Matlab convention
% y(t) + a1 y(t-1) + ... + ap y(t-p) = e(t) + c1 e(t-1) + ... + cq e(t-q):
% AR_1, AR_2, AR_3: the AR coefficients a1, ..., ap (the negatives of the usual AR
%       coefficients)
% MA_1, MA_2: the MA coefficients c1, ..., cq
% maxda, maxdc: the largest estimated standard deviation of the AR and MA
%       coefficients (from the covariance of the parameter estimates)
% noisevar, lossfn, fpe: the noise variance, loss function, and Akaike's final
%       prediction error of the fit
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
% The held-out portion starts at sample floor(pTrain*N), overlapping the training
% portion by one sample.

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
%% Check that a System Identification Toolbox license is available:
% ------------------------------------------------------------------------------
BF_CheckToolbox('identification_toolbox');

% ------------------------------------------------------------------------------
%% Prepare Inputs
% ------------------------------------------------------------------------------
% (1) y, the time series as a column vector
if size(y, 2) > size(y, 1)
	y = y'; % ensure a column vector
end
N = length(y); % number of samples
% Convert y to time series object
y = iddata(y, [], 1);

% orders; vector specifying the AR and MA components
if nargin < 2 || isempty(orders)
	orders = [3, 3]; % AR3, MA3
end
if nargin < 3 || isempty(pTrain)
	pTrain = 0.8; % train on 80% of the data
end
% if nargin < 4 || isempty(trainmode)
%     trainmode = 'first'; % trains on first pTrain proportion of the data.
% end
if nargin < 4 || isempty(numSteps)
	numSteps = 1; % one-step-ahead predictions
end

% ------------------------------------------------------------------------------
%% Fit the model
% ------------------------------------------------------------------------------

% Uses the System Identification Toolbox function armax
m = armax(y, orders);

% ------------------------------------------------------------------------------
%% Statistics on model
% ------------------------------------------------------------------------------

c_ar = m.a; % AR coefficients
c_ma = m.c; % MA coefficients
% Standard deviations of the AR and MA coefficients (the leading, fixed 1 has
% standard deviation 0). These are identical to the dA, dC outputs of
% [A,B,C,D,F,dA,dB,dC] = polydata(m).
da = m.da;
dc = m.dc;

% Make these outputs
if length(c_ar) > 1
	for i = 2:length(c_ar)
		out.(sprintf('AR_%u', i - 1)) = c_ar(i);
	end
end
if length(c_ma) > 1
	for i = 2:length(c_ma)
		out.(sprintf('MA_%u', i - 1)) = c_ma(i);
	end
end

if isempty(da)
	out.maxda = NaN;
else
	out.maxda = max(da);
end
if isempty(dc)
	out.maxdc = NaN;
else
	out.maxdc = max(dc);
end

% ------------------------------------------------------------------------------
% Fit statistics
% ------------------------------------------------------------------------------

% These three measures are basically equivalent -- default hctsa library
% only records fpe.
out.noisevar = m.NoiseVariance; % covariance matrix of noise source
% covmat = m.CovarianceMatrix; % covariance matrix for parameter vector
% parameters = m.ParameterVector; % parameter vector for model: initial values, I'd say...
out.lossfn = m.EstimationInfo.LossFcn;
out.fpe = m.EstimationInfo.FPE; % Final prediction error of model

% out.lastimprovement = m.EstimationInfo.LastImprovement; % Last improvement made in iteration
% (Dropped: aic = aic(m). It is rank-identical to fpe -- Spearman 1.0000 on two
%  collections of real-world series -- so only fpe is kept,
%  matching the choice made in MF_arfit and MF_FitSubsegments.)

% ------------------------------------------------------------------------------
%% Prediction
% ------------------------------------------------------------------------------

% Select first portion of data for estimation
% This could be any portion, actually... Maybe could look at robustness of
% model to different training sets...
ytrain = y(1:floor(pTrain * N));
% ytest = y;
ytest = y(floor(pTrain * N):end); % overlap

% Train the model on just this portion
mp = armax(ytrain, orders);

% Compute step-ahead predictions
% Maybe look at trends across different prediction horizons...
yp = predict(mp, ytest, numSteps, 'init', 'e'); % across whole dataset

mresiduals = yp.y - ytest.y; % prediction minus data (the MF_ResidualAnalysis convention)

% ------------------------------------------------------------------------------
% Get statistics on residuals
% ------------------------------------------------------------------------------
residout = MF_ResidualAnalysis(mresiduals, ytest.y, 'full');

% Convert these to local outputs in quick loop
% Note that default hctsa library does not include rmse field, which is highly
% correlated with the stde field
fields = fieldnames(residout);
for k = 1:length(fields);
	out.(fields{k}) = residout.(fields{k});
end

end
