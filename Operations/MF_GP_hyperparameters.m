% MF_GP_hyperparameters
% 
% Fits a Gaussian Process model to a portion of the time series and returns the
% fitted hyperparameters, as well as statistics describing the goodness of fit
% of the model.
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
% INPUTS:
% y, the input time series
% covfunc, the covariance function, in the standard form fo the gmpl package
% squishorsquash, whether to squash onto the unit interval, or spread across 1:N
% maxN, the maximum length of time series to consider -- greater than this
%               length, time series are resampled down to maxN
% methds, specifies the method of how to resample time series longer than maxN
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

function out = MF_GP_hyperparameters(y,covfunc,squishorsquash,maxN,methds)
% Ben Fulcher, 19/1/2010

%% Preliminaries
doplot = 0; % plot basic outputs
bevocal = 0; % display commentary to command line
N = length(y); % time-series length

%% Check Inputs
if size(y,2) > size(y,1);
    y = y'; % ensure a column vector input
end
% Make sure that y is indeed zscored
if ~BF_iszscored(y)
    warning('The input time series is not, but should be z-scored')
end

if nargin < 2 || isempty(covfunc),
    fprintf(1,'Using a default covariance function: sum of squared exponential and noise\n')
    covfunc = {'covSum', {'covSEiso','covNoise'}};
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
    % specify a proportion of the time series length, N
    maxN = ceil(N*maxN);
end

if nargin < 5 || isempty(methds)
    methds = 'resample';
end

%% Downsample long time series
if (maxN > 0) && (N > maxN)
    switch methds
        case 'resample' % resamples the whole time series down
            f = maxN/N;
            y = resample(y,ceil(f*10000), 10000);
            if length(y) > maxN
                y = y(1:maxN);
            end
            if bevocal
                fprintf(1,'Resampled the time series from a length %u down to %u (%u)\n',N,length(y),maxN);
            end
            N = length(y); % update time series length (should be maxN)
            t = SUB_settimeindex(N,squishorsquash); % set time index

        case 'random_i' % takes maxN random indicies in the time series
            % set time index
            t = SUB_settimeindex(N,squishorsquash);
            % now take samples (unevenly spaced!!)
            ii = randperm(N);
            ii = ii(1:maxN);
            ii = sort(ii,'ascend');
            t = t(ii);
            t = (t-min(t))/max(t)*(maxN-1)+1; % respace from 1:maxN
            y = y(ii);

        case 'random_consec' % takes maxN consecutive indicies from a random position in the time series
            sind = randi(N-maxN+1); % start index
            y = y(sind:sind+maxN-1); % take this bit
            N = length(y); % update time series length (should be maxN)
            t = SUB_settimeindex(maxN,squishorsquash); % set time index
            
        case 'first' % takes first maxN indicies from the time series
            y = y(1:maxN); % take this bit
            N = length(y); % update time series length (should be maxN)
            t = SUB_settimeindex(maxN,squishorsquash); % set time index
            
        case 'random_both' % takes a random starting position and then takes a 1/5 sample from that
            % 1) Take sample from random position in time series
            sind = randi(N-maxN+1); % start index
            y = y(sind:sind+maxN-1); % take this bit
            N = length(y); % update time series length (should be maxN)
            t = SUB_settimeindex(N,squishorsquash); % set time index
            % now take samples (unevenly spaced!!)
            ii = randperm(N);
            ii = ii(1:ceil(maxN/5)); % This 5 is really a parameter...
            ii = sort(ii,'ascend');
            t = t(ii);
            y = y(ii);
        otherwise
            error('Invalid sampling method ''%s''.',methds)
    end
else
    t = SUB_settimeindex(N,squishorsquash); % set time index
end

%% Learn the hyperparameters

% (1) Determine the number of hyperparameters, nhps
s = feval(covfunc{:}); % string in form '2+1', ... tells how many
                        % hyperparameters for each contribution to the
                        % covariance function
nhps = eval(s);

% (2) Intialize hyperparameters before optimization
covfunc1 = covfunc{1};
covfunc2 = covfunc{2};
if strcmp(covfunc1,'covSum') && strcmp(covfunc2{1},'covSEiso') && strcmp(covfunc2{2},'covNoise')
    init_loghyper = ones(3,1)*-1;
    % length parameter is in the ballpark of the difference between time
    % elements
    init_loghyper(1) = log(mean(diff(t)));
else
    init_loghyper = ones(nhps,1)*-1; % Default: initialize all log hyperparameters at -1
end
nfevals = -50; % negative: specifies maximum number of allowed function evaluations


% (3) Perform the optimization
try
    loghyper = minimize(init_loghyper, 'gpr', nfevals, covfunc, t, y);
catch emsg
    if strcmp(emsg.identifier,'MATLAB:posdef')
        disp('Error with lack of positive definite matrix for this function');
        out = NaN; return % return NaN -- the data is not suited to GP fitting
    elseif strcmp(emsg.identifier,'MATLAB:nomem')
        error('Not enough memory to fit a Gaussian Process to this data');
    else
        error('Error fitting Gaussian Process to data, %s\n',emsg.message)
    end
end
hyper = exp(loghyper);

% Output the hyperparameters and log-hyperparameters
for i = 1:nhps
    % Set up structure output
    eval(sprintf('out.h%u = hyper(%u);',i,i));
    eval(sprintf('out.logh%u = loghyper(%u);',i,i));
    % eval(['out.h' num2str(i) ' = hyper(' num2str(i) ');']);
    % eval(['out.logh' num2str(i) ' = loghyper(' num2str(i) ');']);
end


%% For Plotting
if doplot
    xstar = t;
    % xstar = linspace(min(t),max(t),1000)';
    [mu S2] = gpr(loghyper, covfunc, t, y, xstar);
    % S2p = S2 - exp(2*loghyper(3)); % remove noise from predictions
    S2p = S2;
    
    figure('color','w');
    f = [mu+2*sqrt(S2p); flipdim(mu-2*sqrt(S2p),1)];
    fill([xstar; flipdim(xstar,1)], f, [6, 7, 7]/8, 'EdgeColor', [7, 7, 6]/8); % grayscale error bars
    hold on;
    plot(xstar,mu,'k-','LineWidth',2); % mean function
    plot(t,y,'.-k'); % original data
end

%% Other statistics???

% Marginal Likelihood using optimized hyperparameters
out.mlikelihood = - gpr(loghyper, covfunc, t, y);

% Mean error from mean function
[mu, S2] = gpr(loghyper, covfunc, t, y, t); % evaluate at datapoints

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
                            
                            
% Statistics on variance
xstar = linspace(min(t),max(t),1000)'; % crude, I know, but it's nearly 5pm
[mu, S2] = gpr(loghyper, covfunc, t, y, xstar); % evaluate at datapoints
S = sqrt(S2);
out.maxS = max(S); % maximum variance
out.minS = min(S); % minimum variance
out.meanS = mean(S); % mean variance

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

end