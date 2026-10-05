function out = MF_GP_FitAcross(y, covFunc, npoints)
% MF_GP_FitAcross   How well a Gaussian process through a few evenly spaced points reproduces the series.
%
% Trains a Gaussian Process (GP) model, with zero mean and a Gaussian likelihood, on
% points spaced equally throughout the time series, and uses the model to predict all
% the time series values (the intermediate values, and the training points). Times are
% the sample indices. The hyperparameters of the covariance function are learned by
% maximizing the marginal likelihood (MF_GP_LearnHyperp), and the outputs summarize
% the prediction error, the predictive mean and standard deviation, the per-point
% negative log marginal likelihood, and the fitted hyperparameters. If the series is longer than 2000
% samples, predictions are made at 2000 evenly spaced times. A NaN is returned if the
% fit fails.
%
% Uses GP fitting code from the gpml toolbox, which is available here:
% http://gaussianprocess.org/gpml/code.
%
% ---INPUTS:
% y, the input time series
%
% covFunc, the covariance function (structured in the standard way for the gpml
%       toolbox); the default is a sum of a squared-exponential and a noise term,
%       {'covSum',{'covSEiso','covNoise'}}
%
% npoints, the number of points through the time series to fit the GP model to
%       (default 20)
%
% ---OUTPUTS:
% stde, the root-mean-square error of the predictive mean, compared with the series
% meanabs_std, the mean absolute error of the predictive mean, in units of the
%       predictive standard deviation at each time
% stdmu, the standard deviation of the predictive mean over the series
% meanS, stdS, the mean and standard deviation of the predictive standard deviation
%       over the series
% nlml, the negative log marginal likelihood of the whole series (or of the 2000
%       resampled points) under the fitted GP, divided by the number of points, so
%       that it does not grow with the length of the series
% logh1, logh2, logh3, ...: the log hyperparameters of the covariance function, in
%       gpml's order (for the squared-exponential plus noise covariance, the length
%       scale, the signal amplitude, and the noise standard deviation)
% h_lonN, the fitted length scale divided by the series length (only for the
%       squared-exponential plus noise covariance)
%
% ---NOTES:
% In future, the sampling of points could take into account the autocorrelation of
% the time series.

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
%% Check inputs
% ------------------------------------------------------------------------------
if size(y, 2) > size(y, 1);
	y = y'; % make sure a column vector
end

if nargin < 2 || isempty(covFunc)
	fprintf(1, 'Using default sum of SE and noise covariance function\n')
	covFunc = {'covSum', {'covSEiso', 'covNoise'}};
end

if nargin < 3 || isempty(npoints)
	npoints = 20;
end

% ------------------------------------------------------------------------------
%% Get the points
% ------------------------------------------------------------------------------
N = length(y); % time-series length
tt = floor(linspace(1, N, npoints))'; % time range (training)
yt = y(tt);

% ------------------------------------------------------------------------------
%% Optimize the GP parameters for the chosen covariance function
% ------------------------------------------------------------------------------

% Determine the number of hyperparameters, nhps
s = feval(covFunc{:}); % string in form '2+1', ... tells how many
% hyperparameters for each contribution to the
% covariance function
nhps = eval(s);

% Set up the details of the GP:
hyp = struct; % structure for storing hyperparameter information in latest version of GMPL toolbox
meanFunc = {'meanZero'}; hyp.mean = []; % Mean function (mean zero process):
likFunc = @likGauss; hyp.lik = log(0.1); % Likelihood (Gaussian):
% Exact Gaussian inference. (Was @infLaplace: with a Gaussian likelihood the
% Laplace approximation is exact but is computed by Newton iteration, and
% gpml's infLaplace warm-starts that iteration from a PERSISTENT copy of the
% previous call's solution, so every fit depended on whatever series was
% fitted before it -- the same series gave outputs differing in the 4th digit
% from call to call once the optimizer amplified the difference.)
infAlg = @infGaussLik;
nfevals = -50; % number of function evaluations (with negative)

try
	hyp = MF_GP_LearnHyperp(tt, yt, covFunc, meanFunc, likFunc, infAlg, nfevals, hyp);
catch emsg
	error('Error learning hyperparameters for time series')
end
if ~isstruct(hyp) % MF_GP_LearnHyperp returns NaN (not a struct) when the data isn't suited to GP fitting
	out = NaN;
	return
end
loghyper = hyp.cov;
if any(isnan(loghyper))
	out = NaN;
	return
end

% ------------------------------------------------------------------------------
%% Evaluate over the whole space now
% ------------------------------------------------------------------------------
% Evaluate at test points based on training time/data, predicting for
% test times/data
if N <= 2000
	ts = (1:N)';
else % memory constraints force us to crudely resample
	ts = round(linspace(1, N, 2000))';
end
try
	% [mu, S2] = gpr(loghyper, covFunc, tt, yt, ts);
	[mu, S2] = gp(hyp, infAlg, meanFunc, covFunc, likFunc, tt, yt, ts); % evaluate at new time points, ts
catch emsg
	error('Error running Gaussian Process regression on time series: %s', emsg.message);
end

% ------------------------------------------------------------------------------
%% Output statistics
% ------------------------------------------------------------------------------
S = sqrt(S2); % standard deviation function, S
% Root-mean-square error of the mean function, mu.
% (Note: this was previously mean(sqrt((y(ts)-mu).^2)), which cancels pointwise to
%  mean(abs(y(ts)-mu)) -- a mean absolute error, not an RMSE.)
out.stde = sqrt(mean((y(ts) - mu).^2));
out.meanabs_std = mean(abs(y(ts) - mu) ./ S);
out.stdmu = std(mu);
out.meanS = mean(S);
out.stdS = std(S);

% Negative log marginal likelihood per point (gpml's nlZ divided by the number of
% points, so that it does not grow with the number of points, up to 2000)
try
	% out.nlml = - gpr(loghyper, covFunc, ts, y(ts));
	out.nlml = gp(hyp, infAlg, meanFunc, covFunc, likFunc, ts, y(ts)) / length(ts);
catch
	out.nlml = NaN;
end

% Loghyperparameters
for i = 1:nhps
	out.(['logh', num2str(i)]) = loghyper(i); % dynamic field referencing
	% eval(sprintf('out.logh%u = loghyper(%u);',i,i));
end

if strcmp(covFunc{1}, 'covSum') && numel(covFunc{2}) == 2 && ischar(covFunc{2}{1}) && ischar(covFunc{2}{2}) ...
		&& strcmp(covFunc{2}{1}, 'covSEiso') && strcmp(covFunc{2}{2}, 'covNoise')
	% (components with parameters, like {'covMaterniso',3}, are cells, not strings)
	% Give extra output based on length parameter on length of time series
	out.h_lonN = exp(loghyper(1)) / N;
end

%% Subfunctions

%     function loghyper = MF_GP_LearnHyperp(covFunc,nfevals,t,y,init_loghyper)
%         % nfevals--  negative: specifies maximum number of allowed
%         % function evaluations
%         % t: time
%         % y: data
%
%         if nargin < 5 || isempty(init_loghyper)
%             % Use default starting values for parameters
%             % How many hyperparameters
%             s = feval(covFunc{:}); % string in form '2+1', ... tells how many
%             % hyperparameters for each contribution to the
%             % covariance function
%             nhps = eval(s);
%             init_loghyper = ones(nhps,1)*-1; % Initialize all log hyperparameters at -1
%         end
% %         init_loghyper(1) = log(mean(diff(t)));
%
%         % Perform the optimization
%         try
%             loghyper = minimize(init_loghyper, 'gpr', nfevals, covFunc, t, y);
%         catch emsg
%             if strcmp(emsg.identifier,'MATLAB:posdef')
%                 fprintf(1,'Error: lack of positive definite matrix for this function');
%                 loghyper = NaN; return
%             elseif strcmp(emsg.identifier,'MATLAB:nomem')
%                 error('Out of memory');
%                 % return as if a fatal error -- come back to this.
%             else
%                 error('Error fitting Gaussian Process to data')
%             end
%         end
%
%     end

end
