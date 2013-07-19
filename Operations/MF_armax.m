% MF_armax
% 
% Fits an ARMA(p,q) model to the time series and returns various statistics on
% the result.
% 
% Uses the functions iddata, armax, aic, and predict from Matlab's System
% Identification Toolbox
% 
% INPUTS:
% 
% y, the input time series
% 
% orders, a two-vector for p and q, the AR and MA components of the model,
%           respectively,
% 
% ptrain, the proportion of data to train the model on (the remainder is used
%           for testing),
% 
% nsteps, number of steps to predict into the future for testing the model.
% 
% 
% Outputs include the fitted AR and MA coefficients, the goodness of fit in the
% training data, and statistics on the residuals from using the fitted model to
% predict the testing data.

function out = MF_armax(y, orders, ptrain, nsteps)
% Ben Fulcher 1/2/2010

%% Prepare Inputs

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
if nargin < 4 || isempty(nsteps)
    nsteps = 1; % one-step-ahead predictions
end

%% Fit the model
% Uses the system identification toolbox function armax

m = armax(y, orders);

%% Statistics on model

c_ar = m.a; % AR coefficients
c_ma = m.c; % MA coefficients
da = m.da; % must be uncertainties in AR coeffs
dc = m.dc; % must uncertainties in MA coeffs

% Make these outputs
if length(c_ar) > 1
    for i = 2:length(c_ar)
        eval(['out.AR_' num2str(i-1) ' = ' num2str(c_ar(i)) ';']);
    end
end
if length(c_ma) > 1
    for i = 2:length(c_ma)
        eval(['out.MA_' num2str(i-1) ' = ' num2str(c_ma(i)) ';']);
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


% Fit statistics
out.noisevar = m.NoiseVariance; % covariance matrix of noise source
% covmat = m.CovarianceMatrix; % covariance matrix for parameter vector
% parameters = m.ParameterVector; % parameter vector for model: initial values, I'd say...
out.lossfn = m.EstimationInfo.LossFcn;
out.fpe = m.EstimationInfo.FPE; % Final prediction error of model
out.lastimprovement = m.EstimationInfo.LastImprovement; % Last improvement made in interation
out.aic = aic(m); % ~ log(fpe)


%% Prediction

% Select first portion of data for estimation
% This could be any portion, actually... Maybe could look at robustness of
% model to different training sets...
ytrain = y(1:floor(ptrain*N));
% ytest = y;
ytest = y(floor(ptrain*N):end); % overlap

% Train the model on just this portion
mp = armax(ytrain, orders);

% Compute step-ahead predictions
% nsteps = 2; % predicts this many steps ahead
% Maybe look at trends across different prediction horizons...
yp = predict(mp, ytest, nsteps, 'init', 'e'); % across whole dataset

% plot the two:
% plot(y,yp);

mresiduals = ytest.y-yp.y;

% 1) Get statistics on residuals
residout = MF_residanal(mresiduals);

% convert these to local outputs in quick loop
fields = fieldnames(residout);
for k=1:length(fields);
    eval(['out.' fields{k} ' = residout.' fields{k} ';']);
end

% % Train on some proportion, ptrain, of data
% ytrain = y(1:floor(ptrain*N));
% ytest = y(floor(ptrain*N)+1:end);
% mp = armax(ytrain, orders);
% yp = predict(mp,ytrain,2); % 1-step ahead, predictions of length ptrain*N
% % We want only the length of test data worth of predictions, though
% yp = yp(1:N-floor(ptrain*N));
% 
% plot(yp.y,'b');
% hold on;
% plot(ytest.y,'r');

%% This is the way to do it:
% % Select the first half of the data for estimation
% y1 = y(1:48)
% % Estimate a fourth-order autoregressive model
% % using the first half of the data.
% m = ar(y1,4);
% % Compute 6-step ahead prediction
% yhat = predict(m,y,6);
% % Plot the predicted and measured outputs
% plot(y,yhat)


%% Estimate the ARMA model using filter
% e = randn(N,1);
% % yp = filter(c_ma,c_ar,e);
% yp = predict(m,iddata([],e,1),1);
% 
% marma = armax(y,[3,2]);
% compare(y,marma,2)
% 
% %%
% 
% % returns an idpoly model m with estimated parameters and covariances (parameter uncertainties).
% % Estimates the parameters using the prediction-error method and specified orders.
% 
% yp = predict


end