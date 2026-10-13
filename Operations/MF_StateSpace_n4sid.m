function out = MF_StateSpace_n4sid(y, ord, ptrain, steps)
% MF_StateSpace_n4sid   A fitted state-space model of the series, and how well it predicts the later part of the series.
%
% First fits the model to the whole time series, then trains it on the first
% portion and tries to predict the rest.
%
% Model of the form (discrete time, no input, sampling interval 1):
% x(t+1) = A x(t) + K e(t)
% y(t) = C x(t) + e(t)
% for state space matrices A and C, disturbance matrix K (coefficients of
% noise input), vector of states x, and disturbance (noise) e. The model is identified
% by the N4SID subspace algorithm with canonical-variate weighting (BF_SubspaceID).
%
% In the second portion of this code, the state space model is fitted to the
% first p*N samples of the time series, where p is a given proportion and N is
% the length of the time series. This model is then used to predict the latter
% portion of the time series (i.e., the subsequent (1-p)*N samples), and the prediction
% residuals (prediction minus data) are summarized with MF_ResidualAnalysis.
%
% ---INPUTS:
% y, the input time series
% ord, the order of state-space model to implement (can also be the string 'best',
%       for the order from 1 to 10 that minimizes Akaike's information criterion;
%       default: 2)
% ptrain, the proportion of the time series to use for training (default: 0.5)
% steps, the number of steps ahead to predict (default: 1)
%
% ---OUTPUTS:
% From the model fitted to the entire time series:
% A_1, A_2, A_3, A_4, A_5, A_6, A_7, A_8, A_9: the entries of the state-transition
%       matrix A, counted down each column in turn (ord^2 of them)
% k_1, k_2, k_3: the entries of the noise-input vector K (ord of them)
% c_1, c_2, c_3: the entries of the output vector C (ord of them)
% x0mod: the length of the least-squares initial state of the one-step predictor (how
%       far from rest the series starts)
% np: the number of entries of A, K and C (ord^2 + 2*ord)
% Ts: the sampling interval of the model (always 1)
% noisevar: the estimated noise variance, lossfn*N/(N - 2*ord)
% lossfn: the loss function (mean squared one-step prediction error)
% fpe: Akaike's final prediction error, lossfn*(1 + 2*ord/N)/(1 - 2*ord/N)
% bestorder: the order chosen (only when ord = 'best')
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
% ac1diff: absolute lag-1 autocorrelation of the whole time series minus that of
%       the prediction residuals
%
% ---NOTES:
% The entries of A, K, and C are coordinates of the hidden state in the basis of the
% subspace algorithm (the canonical-variate basis), with the sign of each state fixed so
% that C has no negative entries; they are not invariant to a change of state basis.
% The model has 2*ord free parameters once the freedom in the choice of coordinates is
% removed, which is the count used in noisevar and fpe.
% Predictions start from a diffuse (completely unknown) initial state (BF_SSPredict):
% the loss and the residual summaries average over the samples predicted (all but the
% first ord + steps - 1 of the series or of the held-out portion).
% The residuals are prediction minus data (yp - ytest). The held-out portion
% starts at sample floor(ptrain*N), overlapping the training portion by one sample.
% All outputs are NaN if the model cannot be fitted (to the whole series or to the
% training portion), e.g., when the series is too short for the order.

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
%% Check Inputs:
% ------------------------------------------------------------------------------
y = y(:);
N = length(y); % length of time series, N

% (2) the order of the model
if nargin < 2 || isempty(ord)
	ord = 2;
end

% (3) the proportion of the time series to use for training
if nargin < 3 || isempty(ptrain)
	ptrain = 0.5;
end

% (4) the number of steps ahead to predict
if nargin < 4 || isempty(steps)
	steps = 1; % one-step ahead prediction
end

% ------------------------------------------------------------------------------
%% Build the state-space model
% ------------------------------------------------------------------------------
% Use the whole time series -- prediction comes later...
try
	m = BF_SubspaceID(y, ord);
catch
	out = NaN; return % the model cannot be fitted (e.g., too few samples)
end
n = m.order;

if ischar(ord)
	% also return the best order as an output statistic
	out.bestorder = n;
end

% ------------------------------------------------------------------------------
%% Model parameters
% ------------------------------------------------------------------------------
allm_as = m.A(:);
for i = 1:length(allm_as)
	out.(sprintf('A_%u', i)) = allm_as(i);
end
for i = 1:n
	out.(sprintf('k_%u', i)) = m.K(i);
end
for i = 1:n
	out.(sprintf('c_%u', i)) = m.C(i);
end

% Least-squares initial state of the fixed-gain one-step predictor: the predictions are
% affine in the initial state, yhat = yhat0 + B*x0; B(t,i) = C (A - K C)^(t-1) e_i
F = m.A - m.K * m.C;
yhat0 = SUB_FixedGainPredict(m.A, m.K, m.C, y);
B = zeros(N, n);
row = m.C;
for t = 1:N
	B(t, :) = row;
	row = row * F;
end
x0 = BF_Pinv(B) * (y - yhat0);
out.x0mod = sqrt(sum(x0.^2));
out.np = n^2 + 2*n; % the entries of A, K and C

% Transition interval (always 1):
out.Ts = 1;

% Goodness of fit outputs
% Since noisevar, lossfunction, and fpe so highly correlated, the default
% hctsa library only measures fpe.
nEff = 2*n; % free parameters once the choice of state coordinates is removed
out.noisevar = m.loss * N / (N - nEff);
out.lossfn = m.loss;
if nEff >= N
	out.fpe = Inf;
else
	out.fpe = m.loss * (1 + nEff/N) / (1 - nEff/N);
end

% ------------------------------------------------------------------------------
%% Prediction
% ------------------------------------------------------------------------------
% Train the model on the first ptrain proportion, predict the rest (overlapping by one
% sample)
ytrain = y(1:floor(ptrain * N));
ytest = y(floor(ptrain * N):end);
try
	mp = BF_SubspaceID(ytrain, ord);
catch
	out = NaN; return
end
yp = BF_SSPredict(mp.A, mp.K, mp.C, ytest, steps);
mresiduals = yp - ytest; % prediction minus data (the MF_ResidualAnalysis convention)
mresiduals = mresiduals(~isnan(mresiduals)); % the samples that can be predicted

% ------------------------------------------------------------------------------
%% Residual analysis
% ------------------------------------------------------------------------------
residout = MF_ResidualAnalysis(mresiduals, ytest, 'full');

% Convert these to local outputs in quick loop
fields = fieldnames(residout);
for k = 1:length(fields)
	out.(fields{k}) = residout.(fields{k});
end

out.ac1diff = abs(CO_AutoCorr(y, 1, 'Fourier')) - abs(CO_AutoCorr(mresiduals, 1, 'Fourier'));

end
% ------------------------------------------------------------------------------
function yp = SUB_FixedGainPredict(A, K, C, y)
	% One-step predictions of the fixed-gain (steady-state) predictor from the zero state
	N = length(y);
	x = zeros(size(A, 1), 1);
	yp = zeros(N, 1);
	for t = 1:N
		yp(t) = C * x;
		x = A * x + K * (y(t) - yp(t));
	end
end
