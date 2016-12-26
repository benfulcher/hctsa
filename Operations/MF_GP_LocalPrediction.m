function out = MF_GP_LocalPrediction(y,covFunc,numTrain,numTest,numPreds,pmode,randomSeed)
% MF_GP_LocalPrediction     Gaussian Process time-series model for local prediction.
%
% Fits a given Gaussian Process model to a section of the time series and uses
% it to predict to the subsequent datapoint.
%
%---INPUTS:
% y, the input time series
%
% covFunc, covariance function in the standard form for the gpml package.
%           E.g., covFunc = {'covSum', {'covSEiso','covNoise'}} combines squared
%           exponential and noise terms
%
% numTrain, the number of training samples (for each iteration)
%
% numTest, the number of testing samples (for each interation)
%
% numPreds, the number of predictions to make
%
% pmode, the prediction mode:
%       (i) 'beforeafter': predicts the preceding time series values by training
%                           on the following values,
%       (ii) 'frombefore': predicts the following values of the time series by
%                    training on preceding values, and
%       (iii) 'randomgap': predicts random values within a segment of time
%                    series by training on the other values in that segment.
%
% randomSeed, whether (and how) to reset the random seed, using BF_ResetSeed
%               (for 'randomgap' prediction)
%
%---OUTPUTS: summaries of the quality of predictions made, the mean and
% spread of obtained hyperparameter values, and marginal likelihoods.

% Uses GP fitting code from the gpml toolbox, which is available here:
% http://gaussianprocess.org/gpml/code.
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
%% Preliminaries
% ------------------------------------------------------------------------------
doPlot = 0; % plot outputs
N = length(y); % time-series length

% ------------------------------------------------------------------------------
%% Check Inputs
% ------------------------------------------------------------------------------
if size(y,2) > size(y,1)
    y = y'; % ensure a column vector input
end
if nargin < 2 || isempty(covFunc),
    fprintf(1,'Using a default covariance function: sum of squared exponential and noise\n');
    covFunc = {'covSum', {'covSEiso','covNoise'}};
end

if nargin < 3 || isempty(numTrain)
    numTrain = 20; % 20 previous data points to predict the next
end

if nargin < 4 || isempty(numTest)
    numTest = 5; % test on 5 data points into the future
end

if nargin < 5 || isempty(numPreds) % number of predictions
    numPreds = 10; % do it 10 times (equally-spaced) through the time series
end

if nargin < 6 || isempty(pmode)
    pmode = 'frombefore'; % predicts from previous numTrain datapoints
    % can also be 'randomgap' -- fills in random gaps in the middle of a string of
    % data
end

% randomSeed: how to treat the randomization
if nargin < 7
    randomSeed = [];
end

% ------------------------------------------------------------------------------
%% Set up loop
% ------------------------------------------------------------------------------
if ismember(pmode,{'frombefore','randomgap'})
    spns = floor(linspace(1,N-(numTest+numTrain),numPreds)); % starting positions
elseif strcmp(pmode,'beforeafter')
    spns = floor(linspace(1,N-(numTest+numTrain*2),numPreds)); % starting positions
end

% Details of GP:
meanFunc = {'meanZero'}; % zero-mean process
likFunc = @likGauss; % likelihood function (Gaussian)
infAlg = @infLaplace; % Inference algorithm (Laplace approximation)
nfevals = -50;
hyp = struct; % structure for storing hyperparameter information in latest version of GMPL toolbox

% Initialize variables:
mus = zeros(numTest,numPreds); % predicted values
stderrs = zeros(numTest,numPreds); % standard errors on predictions
yss = zeros(numTest,numPreds); % test values
mlikelihoods = zeros(numPreds,1); % marginal likelihoods of model

nhps = eval(feval(covFunc{:})); % number of hyperparameters
loghypers = zeros(nhps,numPreds); % loghyperparameters

for i = 1:numPreds
    %% (0) Set up test and training sets
    switch pmode
    case 'frombefore'
        tt = (1:numTrain)'; % times (make from 1)
        rt = spns(i):spns(i)+numTrain-1; % training range
        yt = y(rt); % training data

        ts = (numTrain+1 : numTrain+1 + numTest-1)'; % times
        rs = spns(i)+numTrain : spns(i)+numTrain + numTest-1; % test range
        ys = y(rs); % test data

    case 'randomgap'
        % Control the random seed (for reproducibility):
        BF_ResetSeed(randomSeed);

        t = (1:numTrain+numTest)';
        r = randperm(numTrain+numTest);
        yy = y(spns(i):spns(i)+numTrain+numTest-1);

        rt = sort(r(1:numTrain),'ascend');
        tt = t(rt);
        yt = yy(rt);

        rs = sort(r(numTrain+1:end),'ascend');
        ts = t(rs);
        ys = yy(rs);

    case 'beforeafter'
        t = (1:2*numTrain+numTest)';
        yy = y(spns(i):spns(i)+2*numTrain+numTest-1);

        rt = [1:numTrain, numTrain+numTest+1:numTrain*2+numTest];
        tt = t(rt);
        yt = yy(rt);

        rs = (numTrain+1 : numTrain+numTest);
        ts = t(rs);
        ys = yy(rs);

    otherwise
        error('Unknown prediction mode ''%s''',pmode);
    end

    % Process to normalize scales
    ys = (ys-mean(yt))/std(yt); % same transformation as training set
    yt = (yt-mean(yt))/std(yt); % zscore training set

    % ------------------------------------------------------------------------------
    %% (1) Learn hyperparameters from training set (t)
    % ------------------------------------------------------------------------------

    % Initialize mean and likelihood
    hyp.mean = []; hyp.lik = log(0.1);
    hyp.cov = [];

    % loghyper = MF_GP_LearnHyperp(covFunc,-50,tt,yt);
    hyp = MF_GP_LearnHyperp(tt,yt,covFunc,meanFunc,likFunc,infAlg,nfevals,hyp);
    loghyper = hyp.cov;

    if isnan(loghyper)
        fprintf(1,'Unable to learn hyperparameters for this time series\n');
        out = NaN; return
    end

    loghypers(:,i) = loghyper;

    % Get marginal likelihood for this model with hyperparameters optimized
    % over training data
    % mlikelihoods(i) = - gpr(loghyper, covFunc, tt, yt);
    mlikelihoods(i) = - gp(hyp, infAlg, meanFunc, covFunc, likFunc, tt, yt);

    % ------------------------------------------------------------------------------
    %% (2) Evaluate at test set (s)
    % ------------------------------------------------------------------------------

    % Evaluate at test points based on training time/data, predicting for
    % test times/data
    % [mu, S2] = gpr(loghyper, covFunc, tt, yt, ts); % old version
    [mu, S2] = gp(hyp, infAlg, meanFunc, covFunc, likFunc, tt, yt, ts); % evaluate at new time points, ts

    % Compare to actual test data --> store in row of errs
    mus(:,i) = mu; % ~predicted values for time series points
    stderrs(:,i) = 2*sqrt(S2); % ~errors on those predictions
    yss(:,i) = ys;

    % Plot
    if doPlot
        if strcmp(pmode,'frombefore')
            plot(tt,yt,'.-k');
            hold on;
            plot(ts,ys,'.-b');
            errorbar(ts,mu,2*sqrt(S2),'m');
            hold off;
        else
            plot(tt,yt,'ok');
            hold on;
            plot(ts,ys,'ob');
            errorbar(ts,mu,2*sqrt(S2),'m');
            hold off;
        end
    end

%     for j=1:numTest
%         % set up structure output
%         err = abs(mu(j)-ys(j))/sqrt(S2(j)); % in units of std at this point
%         eval(['out.abserr' num2str(i) '_' num2str(j) ' = err;']);
%     end

end

% Ok, we're done.

% ------------------------------------------------------------------------------
%% Return statistics on how well it did
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
%% (1) PREDICTION ERROR MEASURES
% ------------------------------------------------------------------------------

% Absolute errors:
allabserrs = abs(mus-yss);
% In units of standard errors (95% confidence interval error bars)
allstderrs = abs(mus-yss)./stderrs;

% ---
% * Stats on all errors:
% ---

% largest error:
out.maxstderr = max(allstderrs(:));
out.maxabserr = max(allabserrs(:));

% smallest error:
out.minstderr = min(allstderrs(:));
out.minabserr = min(allabserrs(:));

% mean error (across all):
out.meanstderr = mean(allstderrs(:));
out.meanabserr = mean(allabserrs(:));

% ---
% * Stats on errors per run
% ---

% Summary of how it did on each run:
stderr_run = mean(allstderrs);
abserr_run = mean(allabserrs);

% Mean error for a run
out.meanstderr_run = mean(stderr_run);
out.meanabserr_run = mean(abserr_run);

% Max error for a run
out.maxstderr_run = max(stderr_run);
out.maxabserr_run = max(abserr_run);

% Min error for a run
out.minstderr_run = min(stderr_run);
out.minabserr_run = min(abserr_run);

% Error bar stats:
out.maxerrbar = max(stderrs(:)); % largest error bar
out.meanerrbar = mean(stderrs(:)); % mean error bar length
out.minerrbar = min(stderrs(:)); % minimum error bar length

% ------------------------------------------------------------------------------
%% (2) HYPERPARAMETER MEASURES
% ------------------------------------------------------------------------------
% mean and std for each hyperparameter
for i = 1:nhps
    out.(sprintf('meanlogh%u',i)) = mean(loghypers(i,:));
    out.(sprintf('stdlogh%u',i)) = std(loghypers(i,:));
end

% ------------------------------------------------------------------------------
%% (3) Marginal likelihood measures
% ------------------------------------------------------------------------------
% Best marginal neg-log-likelihood attained
% Worst marginal neg-log-likelihood attained
% spread in marginal neg-log-likelihoods

out.maxmlik = max(mlikelihoods);
out.minmlik = min(mlikelihoods);
out.stdmlik = std(mlikelihoods);

end
