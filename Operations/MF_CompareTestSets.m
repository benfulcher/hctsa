function out = MF_CompareTestSets(y,theModel,ord,subsetHow,samplep,steps,randomSeed)
% MF_CompareTestSets    Robustness of test-set goodness of fit
%
% Robustness is quantified over different samples in the time series from
% fitting a specified time-series model.
%
% Similar to MF_FitSubsegments, except fits the model on the full time
% series and compares how well it predicts time series in different local
% time-series segments.
%
% Says something of stationarity in spread of values, and something of the
% suitability of model in level of values.
%
% Uses function iddata and predict from Matlab's System Identification Toolbox,
% as well as either ar, n4sid, or armax from Matlab's System Identification
% Toolbox to fit the models, depending on the specified model to fit to the data.
%
%---INPUTS:
% y, the input time series
%
% theModel, the type of time-series model to fit:
%           (i) 'ar', fits an AR model
%           (ii) 'ss', first a state-space model
%           (iii) 'arma', first an ARMA model
%
% ord, the order of the specified model to fit
%
% subsetHow, how to select random subsets of the time series to fit:
%           (i) 'rand', select at random
%           (ii) 'uniform', uniformly distributed segments throughout the time
%                   series
%
% samplep, a two-vector specifying the sampling parameters
%           e.g., [20, 0.1] repeats 20 times for segments 10% the length of the
%                           time series
%
% steps, the number of steps ahead to do the predictions.
%
% randomSeed, whether (and how) to reset the random seed, using BF_ResetSeed
%               (when 'rand' specified for subsetHow)

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
%% Check that a System Identification Toolbox license is available:
% ------------------------------------------------------------------------------
BF_CheckToolbox('identification_toolbox')

% ------------------------------------------------------------------------------
%% Preliminaries
% ------------------------------------------------------------------------------
N = length(y); % length of time series

% ------------------------------------------------------------------------------
%% Check inputs, set defaults
% ------------------------------------------------------------------------------
% (1) y: column vector time series
if nargin < 1 || isempty(y)
    error('No input time series provided');
end

% Convert y to time series object for the System Identification Toolbox
y = iddata(y,[],1);

% (2) Model, the type of model to fit
if nargin < 2 || isempty(theModel)
    theModel = 'ss';
    % Fit a state space model by default
end

% (3) Model order, ord
if nargin < 3 || isempty(ord)
    ord = 2;
    % model of order 2 by default. Not the best defaults.
end

% (4) How to choose subsets from the time series, subsetHow
if nargin < 4 || isempty(subsetHow)
    subsetHow = 'rand'; % takes segments randomly from time series
end

% (5) Sampling parameters, samplep
if nargin < 5 || isempty(samplep)
    samplep = [20, 0.1]; % sample 20 times with 10%-length subsegments
end

% (6) Predict some number of steps ahead in test sets, steps
if nargin < 6 || isempty(steps)
    steps = 2; % default: predict 2 steps ahead in test set
end

% (7)  randomSeed: how to treat the randomization
if nargin < 7
    randomSeed = [];
end

% ------------------------------------------------------------------------------
%% Fit the model
% ------------------------------------------------------------------------------
% model will be stored as a model object, m
% model is fitted using the entire dataset as the training set
% test sets will be smaller chunks of this.
% [could also fit multiple models using data not in multiple test sets, but
% this is messier]

switch theModel
    case 'ar' % fit an ar model of specified order
        if strcmp(ord,'best')
            % Use arfit software to retrieve the optimum ar order by some
            % criterion (Schwartz's Bayesian Criterion, SBC)
            try
                [west, Aest, Cest, SBC, FPE, th] = ARFIT_arfit(y.y, 1, 10, 'sbc', 'zero');
            catch
                error('Error running ''arfit'' -- is the ARFIT toolbox installed?')
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
        error('Unknown model ''%s''', theModel);
end

% ------------------------------------------------------------------------------
%% Prepare to do a series of predictions
% ------------------------------------------------------------------------------
% Number of samples to take, numPred
numPred = samplep(1);
% Initialize quantities to store into
rmserrs = zeros(numPred,1);
mabserrs = zeros(numPred,1);
ac1s = zeros(numPred,1);
meandiffs = zeros(numPred,1);
stdrats = zeros(numPred,1);

% Set ranges beforehand
r = zeros(numPred,2);

switch subsetHow
    case 'rand'
        if samplep(2) < 1 % specified a fraction of time series
            l = floor(N*samplep(2));
        else % specified an absolute interval
            l = samplep(2);
        end

        % Control the random seed (for reproducibility):
        BF_ResetSeed(randomSeed);

        % numPred starting points:
        spts = randi(N-l+1,numPred,1);
        r(:,1) = spts;
        r(:,2) = spts+l-1;

    case 'uniform'
        if length(samplep) == 1 % size will depend on number of unique subsegments
            spts = round(linspace(0,N,numPred+1)); % numPred+1 boundaries = numPred portions
            r(:,1) = spts(1:numPred)+1;
            r(:,2) = spts(2:end);
        else
            if samplep(2) < 1 % specified a fraction of time series
                l = floor(N*samplep(2));
            else % specified an absolute interval
                l = samplep(2);
            end
            spts = round(linspace(1,N-l+1,numPred)); % numPred+1 boundaries = numPred portions
            r(:,1) = spts;
            r(:,2) = spts+l-1;
        end
    otherwise
        error('Unknown subset method ''%s''',subsetHow);
end

% Quickly check that ranges are valid
if any(r(:,1) >= r(:,2))
    error('Invalid settings');
end

% ------------------------------------------------------------------------------
%% Do the series of predictions
% ------------------------------------------------------------------------------
for i = 1:numPred
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
    ac1s(i) = CO_AutoCorr(mres,1,'Fourier');

    % Get statistics on output time series
    meandiffs(i) = abs(mean(yp.y) - mean(ytest.y));
    stdrats(i) = abs(std(yp.y)/std(ytest.y));

%     % 1) Get statistics on residuals
%     residout = MF_ResidualAnalysis(mresiduals);
%
%     % convert these to local outputs in quick loop
%     fields = fieldnames(residout);
%     for k=1:length(fields);
%         eval(['out.' fields{k} ' = residout.' fields{k} ';']);
%     end
end

% ------------------------------------------------------------------------------
%% Return statistics on outputs
% ------------------------------------------------------------------------------

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
