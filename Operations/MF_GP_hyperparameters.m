function out = MF_GP_hyperparameters(y,covFunc,squishorsquash,maxN,resampleHow,randomSeed)
% MF_GP_hyperparameters    Gaussian Process time-series model parameters and goodness of fit
%
% Uses GP fitting code from the gpml toolbox, which is available here:
% http://gaussianprocess.org/gpml/code.
%
% The code can accomodate a range of covariance functions, e.g.:
% (i) a sum of squared exponential and noise terms, and
% (ii) a sum of squared exponential, periodic, and noise terms.
%
% The model is fitted to <> samples from the time series, which are
% chosen by:
% (i) resampling the time series down to this many data points,
% (ii) taking the first 200 samples from the time series, or
% (iii) taking random samples from the time series.
%
%---INPUTS:
% y, the input time series
%
% covFunc, the covariance function, in the standard form of the gmpl package
%
% squishorsquash, whether to squash onto the unit interval, or spread across 1:N
%
% maxN, the maximum length of time series to consider -- inputs greater than
%           this length are resampled down to maxN
%
% resampleHow, specifies the method of how to resample time series longer than maxN
%
% randomSeed, whether (and how) to reset the random seed, using BF_ResetSeed,
%             for settings of resampleHow that involve random number generation

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
doPlot = 0; % plot basic outputs
beVocal = 0; % display commentary to command line
N = length(y); % time-series length

% ------------------------------------------------------------------------------
%% Check Inputs
% ------------------------------------------------------------------------------
if size(y,2) > size(y,1);
    y = y'; % ensure a column vector input
end
% Make sure that y is indeed zscored
if ~BF_iszscored(y)
    warning('The input time series is not, but should be z-scored')
end

if nargin < 2 || isempty(covFunc),
    fprintf(1,'Using a default covariance function: sum of squared exponential and noise\n');
    covFunc = {'covSum', {'covSEiso','covNoise'}};
end

if nargin < 3 || isempty(squishorsquash)
    squishorsquash = 1;
end

if nargin < 4 || isempty(maxN)
    maxN = 500; % maximum length time series we do this for --
                 % resample longer time series
    % maxN = 0 --> include the whole thing
end
if (maxN > 0) && (maxN < 1)
    % Specify a proportion of the time series length, N
    maxN = ceil(N*maxN);
end

if nargin < 5 || isempty(resampleHow)
    resampleHow = 'resample';
end

if nargin < 6
    randomSeed = [];
end

% Inference algorithm -- use the Laplace approximation:
infAlg = @infLaplace;

% ------------------------------------------------------------------------------
%% Downsample long time series
% ------------------------------------------------------------------------------
if (maxN > 0) && (N > maxN)
    switch resampleHow
        case 'resample' % resamples the whole time series down
            f = maxN/N;
            y = resample(y,ceil(f*10000), 10000);
            if length(y) > maxN
                y = y(1:maxN);
            end
            if beVocal
                fprintf(1,'Resampled the time series from a length %u down to %u (%u)\n',N,length(y),maxN);
            end
            N = length(y); % update time series length (should be maxN)
            t = SUB_settimeindex(N,squishorsquash); % set time index

        case 'random_i' % takes maxN random indicies in the time series
            % Set time index
            t = SUB_settimeindex(N,squishorsquash);
            % Control the random seed (for reproducibility):
            BF_ResetSeed(randomSeed);
            % Now take samples (unevenly spaced!!)
            ii = randsample(N,maxN);
            ii = sort(ii,'ascend');
            t = t(ii);
            t = (t-min(t))/max(t)*(maxN-1)+1; % respace from 1:maxN
            y = y(ii);

        case 'random_consec' % takes maxN consecutive indicies from a random position in the time series
            % Control the random seed (for reproducibility):
            BF_ResetSeed(randomSeed);
            sind = randi(N-maxN+1); % start index
            y = y(sind:sind+maxN-1); % take this bit
            t = SUB_settimeindex(maxN,squishorsquash); % set time index

        case 'first' % takes first maxN indicies from the time series
            y = y(1:maxN); % take this bit
            t = SUB_settimeindex(maxN,squishorsquash); % set time index

        case 'random_both' % takes a random starting position and then takes a 1/5 sample from that
            % Control the random seed (for reproducibility):
            BF_ResetSeed(randomSeed);
            % Take sample from random position in time series
            sind = randi(N-maxN+1); % start index
            y = y(sind:sind+maxN-1); % take this bit
            N = length(y); % update time series length (should be maxN)
            t = SUB_settimeindex(N,squishorsquash); % set time index
            % Now take samples (unevenly spaced!!)
            ii = randsample(N,ceil(maxN/5)); % This 5 is really a parameter...
            ii = sort(ii,'ascend');
            t = t(ii);
            y = y(ii);

        otherwise
            error('Invalid sampling method ''%s''.',resampleHow)
    end
else
    t = SUB_settimeindex(N,squishorsquash); % set time index
end

% ------------------------------------------------------------------------------
%% Learn the hyperparameters
% ------------------------------------------------------------------------------

% (1) Determine the number of hyperparameters, numHPs
s = feval(covFunc{:}); % string in form '2+1', ... tells how many
                        % hyperparameters for each contribution to the
                        % covariance function
numHPs = eval(s); % number of hyperparameters

% (2) Intialize hyperparameters before optimization and perform the optimization
hyp = struct; % structure for storing hyperparameter information in latest version of GMPL toolbox

% Mean function (mean zero process):
meanFunc = {'meanZero'}; hyp.mean = [];

% Likelihood (Gaussian):
likFunc = @likGauss; hyp.lik = log(0.1);

% Maximum number of allowed function evaluations
numfevals = -50; % (specified as the negative)

hyp = MF_GP_LearnHyperp(t,y,covFunc,meanFunc,likFunc,infAlg,numfevals,hyp);

% Get non-logarithmic hyperparameters
logHyper = hyp.cov;
hyper = exp(logHyper);

% Output the hyperparameters and log-hyperparameters
for i = 1:numHPs
    % Set up structure output
    out.(sprintf('h%u',i)) = hyper(i);
    out.(sprintf('logh%u',i)) = logHyper(i);
end

% ------------------------------------------------------------------------------
%% For Plotting
% ------------------------------------------------------------------------------
if doPlot
    xstar = t;
    % xstar = linspace(min(t),max(t),1000)';
    [mu, S2] = gpr(logHyper, covFunc, t, y, xstar);
    % S2p = S2 - exp(2*logHyper(3)); % remove noise from predictions
    S2p = S2;

    figure('color','w');
    f = [mu+2*sqrt(S2p); flipdim(mu-2*sqrt(S2p),1)];
    fill([xstar; flipdim(xstar,1)], f, [6, 7, 7]/8, 'EdgeColor', [7, 7, 6]/8); % grayscale error bars
    hold on;
    plot(xstar,mu,'k-','LineWidth',2); % mean function
    plot(t,y,'.-k'); % original data
end

% ------------------------------------------------------------------------------
%% Other statistics???
% ------------------------------------------------------------------------------

% Negative log marginal likelihood using optimized hyperparameters
out.mlikelihood = gp(hyp, infAlg, meanFunc, covFunc, likFunc, t, y);

% Mean error from fit
[mu, S2] = gp(hyp, infAlg, meanFunc, covFunc, likFunc, t, y, t); % evaluate at datapoints
% [mu, S2] = gpr(logHyper, covFunc, t, y, t); % evaluate at datapoints

if std(mu) < 0.01; % hasn't fit the time series well at all -- too constant
    fprintf(1,'This time series is not suited to Gaussian Process fitting\n');
    out = NaN; return
end

out.rmserr = mean(sqrt((y-mu).^2));
% Better to look at mean distance away in units of std
out.mabserr_std = mean(abs((y-mu)./sqrt(S2)));
out.std_mu_data = std(mu); % std of mean function evaluated at datapoints
                            % (if not close to one, means a problem with
                            % fitting)
out.std_S_data = std(sqrt(S2)); % should vary a fair bit


% Statistics on variance:
xstar = linspace(min(t),max(t),1000)'; % crude, I know, but it's nearly 5pm
[~, S2] = gpr(logHyper, covFunc, t, y, xstar); % evaluate at datapoints
S = sqrt(S2);
out.maxS = max(S); % maximum variance
out.minS = min(S); % minimum variance
out.meanS = mean(S); % mean variance

% ------------------------------------------------------------------------------
function t = SUB_settimeindex(N,squishorsquash)
    %% Set time index
    % Difficult for processes on different time scales -- to squash them all
    % into one time 'window' with linspace, or spread them all out into a
    % single 'sampling rate' with 1:N...?
    if squishorsquash
        t = (1:N)';
    else
        t = linspace(0,1,N)';
    end
end
% ------------------------------------------------------------------------------

end
