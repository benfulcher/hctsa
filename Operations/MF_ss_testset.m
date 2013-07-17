function out = MF_ss_testset(y,model,ord,howtosubset,samplep,steps)
% Looks at robustness of test set goodness of fit over different samples in
% the time series from fitting a model (of given order or
% optimal: 'best'). Something of stationarity in spread of values.
% Something of suitability of model in level of values.
% Uses function iddata and predict from Matlab's System Identification Toolbox
% Also uses either ar, n4sid, or armax from Matlab's System Identification Toolbox to fit the models, depending on the model specified as the input, model
% Ben Fulcher 12/2/2010

%% Preliminaries
N = length(y); % length of time series

%% Inputs

% (1) y: column vector time series
if nargin < 1 || isempty(y)
    error('No input time series provided');
end

% Convert y to time series object for the System Identification Toolbox
y = iddata(y,[],1);

% (2) Model, the type of model to fit
if nargin < 2 || isempty(model)
    model = 'ss';
    % Fit a state space model by default
end
    
% (3) Model order, ord
if nargin < 3 || isempty(ord)
    ord = 2;
    % model of order 2 by default. Not the best defaults.
end

% (4) How to choose subsets from the time series, howtosubset
if nargin < 4 || isempty(howtosubset)
    howtosubset = 'rand'; % takes segments randomly from time series
end

% (5) Sampling parameters, samplep
if nargin < 5 || isempty(samplep)
    samplep = [20, 0.1]; % sample 20 times with 10%-length subsegments
end

% (6) Predict some number of steps ahead in test sets, steps
if nargin < 6 || isempty(steps)
    steps = 2; % default: predict 2 steps ahead in test set
end

%% Fit the model
% model will be stored as a model object, m
% model is fitted using the entire dataset as the training set
% test sets will be smaller chunks of this.
% [could also fit multiple models using data not in multiple test sets, but
% this is messier]
switch model
    case 'ar' % fit an ar model of specified order
        if strcmp(ord,'best')
            % Use arfit software to retrieve the optimum ar order by some
            % criterion (Schwartz's Bayesian Criterion, SBC)
            try
                [west, Aest, Cest, SBC, FPE, th] = arfit(y.y, 1, 10, 'sbc', 'zero');
            catch
                error('Error running ''arfit'' -- have you installed the ARFIT toolbox?')
            end
            ord = length(Aest);
        end
        m = ar(y,ord);
    case 'ss' % fit a state space model of specified order
        m = n4sid(y,ord);
    case 'arma' % fit an arma model of specified orders
        % Note: order should be a two-component vector
        m = armax(y,ord);
    otherwise
        error('Unknown model ''%s''', model);
end

%% Prepare to do a series of predictions
% Number of samples to take, npred
npred = samplep(1);
% Initialize quantities to store into
rmserrs = zeros(npred,1);
mabserrs = zeros(npred,1);
ac1s = zeros(npred,1);
meandiffs = zeros(npred,1);
stdrats = zeros(npred,1);

% Set ranges beforehand
r = zeros(npred,2);
switch howtosubset
    case 'rand'
        if samplep(2) < 1 % specified a fraction of time series
            l = floor(N*samplep(2));
        else % specified an absolute interval
            l = samplep(2);
        end
        spts = randi(N-l+1,npred,1); % npred starting points
        r(:,1) = spts;
        r(:,2) = spts+l-1;
        
    case 'uniform'
        if length(samplep) == 1 % size will depend on number of unique subsegments
            spts = round(linspace(0,N,npred+1)); % npred+1 boundaries = npred portions
            r(:,1) = spts(1:npred)+1;
            r(:,2) = spts(2:end);
        else
            if samplep(2)<1 % specified a fraction of time series
                l = floor(N*samplep(2));
            else % specified an absolute interval
                l = samplep(2);
            end
            spts = round(linspace(1,N-l+1,npred)); % npred+1 boundaries = npred portions
            r(:,1) = spts;
            r(:,2) = spts+l-1;
        end
    otherwise
        error('Unknown subset method ''%s''',howtosubset);
end

% Quickly check that ranges are valid
if any(r(:,1) >= r(:,2))
    error('Invalid settings');
end

%% Do the series of predictions
for i = 1:npred
    % retrieve the test set

    ytest = y(r(i,1):r(i,2));
    
    % Compute step-ahead predictions using System Identification Toolbox:
    yp = predict(m, ytest, steps); % across test set using model, m,
                                   % fitted to entire data set

%     e = pe(m, ytest); % prediction errors -- exactly the same as returning 
%                       % residuals of 1-step-ahead prediction model


    % plot the two:
%     plot(y,yp);

    % Get statistics on residuals
    mres = yp.y - ytest.y;
    
    rmserrs(i) = sqrt(mean(mres.^2));
    mabserrs(i) = mean(abs(mres));
    ac1s(i) = CO_autocorr(mres,1);
    
    % Get statistics on output time series
    meandiffs(i) = abs(mean(yp.y) - mean(ytest.y));
    stdrats(i) = abs(std(yp.y))/abs(std(ytest.y));
    
    
%     % 1) Get statistics on residuals
%     residout = MF_residanal(mresiduals);
% 
%     % convert these to local outputs in quick loop
%     fields = fieldnames(residout);
%     for k=1:length(fields);
%         eval(['out.' fields{k} ' = residout.' fields{k} ';']);
%     end

    
end

%% Return statistics on outputs

% RMS errors, rmserrs
out.rmserr_mean = mean(rmserrs);
out.rmserr_median = median(rmserrs);
out.rmserr_std = std(rmserrs);
out.rmserr_iqr = iqr(rmserrs);

% mean absolute errors, mabserrs
out.mabserr_mean = mean(mabserrs);
out.mabserr_median = median(mabserrs);
out.mabserr_std = std(mabserrs);
out.mabserr_iqr = iqr(mabserrs);

% Autocorrelations at lag 1, ac1s
% NOT absolute values of ac1s... absolute values of operations on *raw*
% ac1s...
out.ac1s_mean = abs(mean(ac1s));
out.ac1s_median = abs(median(ac1s));
out.ac1s_std = std(ac1s);   
out.ac1s_iqr = iqr(ac1s);

% Differences in mean between two series
out.meandiffs_mean = mean(meandiffs);
out.meandiffs_median = median(meandiffs);
out.meandiffs_std = std(meandiffs);
out.meandiffs_iqr = iqr(meandiffs);

% Ratio of standard deviations between two series
out.stdrats_mean = mean(stdrats);
out.stdrats_median = median(stdrats);
out.stdrats_std = std(stdrats);
out.stdrats_iqr = iqr(stdrats);


end