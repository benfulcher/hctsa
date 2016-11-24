function out = MF_armax(y, orders, ptrain, numSteps)
% MF_armax  Statistics on a fitted ARMA model.
%
% Uses the functions iddata, armax, aic, and predict from Matlab's System
% Identification Toolbox
%
%---INPUTS:
%
% y, the input time series
%
% orders, a two-vector for p and q, the AR and MA components of the model,
%           respectively,
%
% ptrain, the proportion of data to train the model on (the remainder is used
%           for testing),
%
% numSteps, number of steps to predict into the future for testing the model.
%
%
%---OUTPUTS: include the fitted AR and MA coefficients, the goodness of fit in
% the training data, and statistics on the residuals from using the fitted model
% to predict the testing data.

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
%% Check that a System Identification Toolbox license is available:
% ------------------------------------------------------------------------------
BF_CheckToolbox('identification_toolbox')

% ------------------------------------------------------------------------------
%% Prepare Inputs
% ------------------------------------------------------------------------------
% (1) y, the time series as a column vector
if size(y,2) > size(y,1)
   y = y'; % ensure a column vector
end
N = length(y); % number of samples
% Convert y to time series object
y = iddata(y,[],1);

% orders; vector specifying the AR and MA components
if nargin < 2 || isempty(orders)
    orders = [3, 3]; % AR3, MA3
end
if nargin < 3 || isempty(ptrain)
    ptrain = 0.8; % train on 80% of the data
end
% if nargin < 4 || isempty(trainmode)
%     trainmode = 'first'; % trains on first ptrain proportion of the data.
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
da = m.da; % must be uncertainties in AR coeffs
dc = m.dc; % must uncertainties in MA coeffs

% Make these outputs
if length(c_ar) > 1
    for i = 2:length(c_ar)
        out.(sprintf('AR_%u',i-1)) = c_ar(i);
    end
end
if length(c_ma) > 1
    for i = 2:length(c_ma)
        out.(sprintf('MA_%u',i-1)) = c_ma(i);
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

out.lastimprovement = m.EstimationInfo.LastImprovement; % Last improvement made in interation
out.aic = aic(m); % ~ log(fpe)

% ------------------------------------------------------------------------------
%% Prediction
% ------------------------------------------------------------------------------

% Select first portion of data for estimation
% This could be any portion, actually... Maybe could look at robustness of
% model to different training sets...
ytrain = y(1:floor(ptrain*N));
% ytest = y;
ytest = y(floor(ptrain*N):end); % overlap

% Train the model on just this portion
mp = armax(ytrain, orders);

% Compute step-ahead predictions
% Maybe look at trends across different prediction horizons...
yp = predict(mp, ytest, numSteps, 'init', 'e'); % across whole dataset

mresiduals = ytest.y - yp.y;

% ------------------------------------------------------------------------------
% Get statistics on residuals
% ------------------------------------------------------------------------------
residout = MF_ResidualAnalysis(mresiduals);

% Convert these to local outputs in quick loop
% Note that default hctsa library does not include rmse field, which is highly
% correlated with the stde field
fields = fieldnames(residout);
for k = 1:length(fields);
    out.(fields{k}) = residout.(fields{k});
end


end
