function out = MF_sissm(y,order,ptrain,steps)
% fits a state-space model using methods in system identification toolbox
% of MATLAB
% First fits model to the whole time series; then trains it on the first
% portion and tries to predict the rest.
% Ben Fulcher 1/2/2010

%% Inputs

% (1) y: the time series as a column vector
% Convert y to time series object
N = length(y); % length of time series, N
y = iddata(y,[],1);

% (2) Order, the order of the state space model to fit. Can specify a positive
% integer or the string 'best'.
if nargin<2 || isempty(order)
    order = 2;
end

% (3) train model on this proportion.
if nargin<3 || isempty(ptrain)
    ptrain = 0.5;
end

% (4) steps, step-ahead prediction
if nargin<4 || isempty(steps)
    steps = 1; % one-step ahead prediction
end




%% Build the state-space model
% use the whole time series -- prediction comes later...
% m = n4sid(y,'best'); % chooses the best model order among orders 1:10
m = n4sid(y,order); % fits a model of given order

if strcmp(order,'best')
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
for i=1:length(allm_as)
    eval(['out.A_' num2str(i) ' = ' num2str(allm_as(i)) ';']);
end
for i=1:length(m_ks)
    eval(['out.k_' num2str(i) ' = ' num2str(m_ks(i)) ';']);
end
for i=1:length(m_cs)
    eval(['out.c_' num2str(i) ' = ' num2str(m_cs(i)) ';']);
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
% keyboard


%% Prediction

% Select first portion of data for estimation
% This could be any portion, actually... Maybe could look at robustness of
% model to different training sets...
ytrain = y(1:floor(ptrain*N));
% ytest = y;
ytest = y(floor(ptrain*N):end); % overlap

% Train the model on just this portion
% mp = armax(ytrain, orders);
try mp = n4sid(ytrain, order);
catch emsg
    disp('couldn''t fit the model, I''m afraid')
    out = NaN; return
end

% Compute step-ahead predictions
% steps = 2; % predicts this many steps ahead
% Maybe look at trends across different prediction horizons...
yp = predict(mp, ytest, steps, 'init', 'e'); % across whole ytest dataset

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

out.ac1diff = abs(CO_autocorr(y.y,1))-abs(CO_autocorr(mresiduals,1));

end