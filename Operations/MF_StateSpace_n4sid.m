% MF_StateSpace_n4sid
% 
% Fits a state space model of given order to the time series using the n4sid
% function in Matlab's System Identification Toolbox.
% 
% First fits the model to the whole time series, then trains it on the first
% portion and tries to predict the rest.
% 
% In the second portion of this code, the state space model is fitted to the
% first p*N samples of the time series, where p is a given proportion and N is
% the length of the time series.
% 
% This model is then used to predict the latter portion of the time
% series (i.e., the subsequent (1-p)*N samples).
% 
% Uses the functions iddata, n4sid, aic, and predict from Matlab's System Identification Toolbox
% 
% INPUTS:
% y, the input time series
% ord, the order of state-space model to implement (can also be the string 'best')
% ptrain, the proportion of the time series to use for training
% steps, the number of steps ahead to predict
% 
% Outputs are parameters from the model fitted to the entire time series, and
% goodness of fit and residual analysis from n4sid prediction.
% 
% ------------------------------------------------------------------------------
% Copyright (C) 2013,  Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones., "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2010). DOI: 10.1098/rsif.2013.0048
%
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
% details.
% 
% You should have received a copy of the GNU General Public License along with
% this program.  If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

function out = MF_StateSpace_n4sid(y,ord,ptrain,steps)
% Ben Fulcher, 1/2/2010

%% Check Inputs:

% (1) y: the time series as a column vector
% Convert y to time series object
N = length(y); % length of time series, N

if ~exist('iddata')
    error('''iddata'' not found, you''ll need the System Identification Toolbox to run this code');
end
y = iddata(y,[],1);

% (2) Order, the order of the state space model to fit. Can specify a positive
% integer or the string 'best'.
if nargin < 2 || isempty(ord)
    ord = 2;
end

% (3) train model on this proportion.
if nargin < 3 || isempty(ptrain)
    ptrain = 0.5;
end

% (4) steps, step-ahead prediction
if nargin < 4 || isempty(steps)
    steps = 1; % one-step ahead prediction
end


%% Build the state-space model
% use the whole time series -- prediction comes later...
% m = n4sid(y,'best'); % chooses the best model order among orders 1:10
m = n4sid(y,ord); % fits a state-space model of given order

if strcmp(ord,'best')
    % also return the best order as an output statistic
    out.bestorder = length(m.k);
end

%% Model parameters
% Analysis of model

% Parameters
m_as = m.A; % 'transition' matrix in underlying ss model, x
m_ks = m.K; % coefficients for noise input in ss model
m_cs = m.C; % coefficients for measurement function, y
m_x0 = m.X0; % initial condition
m_np = length(m.ParameterVector); % number of parameters fitted

% output parameters
allm_as = m_as(:);
for i = 1:length(allm_as)
    eval(sprintf('out.A_%u = allm_as(%u);',i,i));
end
for i = 1:length(m_ks)
    eval(sprintf('out.k_%u = m_ks(%u);',i,i));
end
for i = 1:length(m_cs)
    eval(sprintf('out.c_%u = m_cs(%u);',i,i));
end
out.x0mod = sqrt(sum(m_x0.^2));
out.np = m_np;


% Goodness of fit outputs
out.m_noisevar = m.NoiseVariance; % a scalar number
out.m_Ts = m.Ts; % Transition interval
out.m_lossfn = m.EstimationInfo.LossFcn;
out.m_fpe = m.EstimationInfo.FPE;
out.m_aic = aic(m);
% out.m_bic = out.m_aic - 2*m_np + m_np*log(N); % scrappy and possibly wrong

% plot(y)

%% Prediction

% Select first portion of data for estimation
% This could be any portion, actually... Maybe could look at robustness of
% model to different training sets...
ytrain = y(1:floor(ptrain*N));
% ytest = y;
ytest = y(floor(ptrain*N):end); % overlap

% Train the model on just this portion
% mp = armax(ytrain, orders);
try
    mp = n4sid(ytrain, ord);
catch emsg
    error('Couldn''t fit the model to this time series: %s',emsg.message)
    % out = NaN; return
end

% Compute step-ahead predictions
% steps = 2; % predicts this many steps ahead
% Maybe look at trends across different prediction horizons...
yp = predict(mp, ytest, steps, 'init', 'e'); % across whole ytest dataset

% plot the two:
% plot(y,yp);

mresiduals = ytest.y-yp.y;

% 1) Get statistics on residuals
residout = MF_ResidualAnalysis(mresiduals);

% convert these to local outputs in quick loop
fields = fieldnames(residout);
for k = 1:length(fields);
    eval(sprintf('out.%s = residout.%s;',fields{k},fields{k}));
end

out.ac1diff = abs(CO_AutoCorr(y.y,1))-abs(CO_AutoCorr(mresiduals,1));

end