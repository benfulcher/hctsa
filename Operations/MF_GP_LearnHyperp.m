function hyp = MF_GP_LearnHyperp(t,y,covFunc,meanFunc,likFunc,infAlg,nfevals,hyp)
% MF_GP_LearnHyperp     Learns Gaussian Process hyperparameters for a time series
%
% Used by main Gaussian Process model fitting operations.
%
% References code 'minimize' from the GAUSSIAN PROCESS REGRESSION AND
% CLASSIFICATION Toolbox version 3.2, which is avilable at:
% http://gaussianprocess.org/gpml/code
%
%---INPUTS:
%
% t,             time
% y,             data
% covFunc,       the covariance function, formatted as gpml likes it
% meanFunc, the mean function, formatted as gpml likes it
% likFunc, the likelihood function, formatted as gpml likes it
% infAlg, the inference algorithm (in gpml form)
% nfevals,       the number of function evaluations

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

if nargin < 3 || isempty(covFunc)
    covFunc = @covSEiso;
end
if nargin < 4 || isempty(meanFunc)
    % Mean function (mean zero process):
    meanFunc = {'meanZero'}; hyp.mean = [];
end
if nargin < 5 || isempty(likFunc)
    likFunc = @likGauss; % negative: specifies maximum number of allowed function evaluations
    hyp.lik = log(0.1);
end
if nargin < 6 || isempty(infAlg)
    % Inference algorithm:
    infAlg = @infLaplace;
end
if nargin < 7 || isempty(nfevals)
    nfevals = -50; % negative: specifies maximum number of allowed function evaluations
end
% ------------------------------------------------------------------------------

% Number of hyperparameters:
s = feval(covFunc{:});
nhps = eval(s);

% Initial values for covariance function:
covFunc1 = covFunc{1};
covFunc2 = covFunc{2};
if strcmp(covFunc1,'covSum') && strcmp(covFunc2{1},'covSEiso') && strcmp(covFunc2{2},'covNoise')
    hyp.cov = zeros(3,1);
    % length parameter is in the ballpark of the difference between time
    % elements
    hyp.cov(1) = log(mean(diff(t)));
else
    hyp.cov = zeros(nhps,1); % Default: initialize all log hyperparameters at -1
end

% ------------------------------------------------------------------------------
% Perform the optimization
% ------------------------------------------------------------------------------
try
    % loghyper = minimize(init_loghyper, 'gpr', nfevals, covFunc, t, y);
    hyp = minimize(hyp, @gp, nfevals, infAlg, meanFunc, covFunc, likFunc, t, y);
catch emsg
    if strcmp(emsg.identifier,'MATLAB:posdef')
        fprintf(1,'Error with lack of positive definite matrix for this function\n');
        hyp = NaN; return % return NaN -- the data is not suited to GP fitting
    elseif strcmp(emsg.identifier,'MATLAB:nomem')
        error('Not enough memory to fit a Gaussian Process to this data');
    else
        error('Error fitting Gaussian Process to data: %s\n',emsg.message)
    end
end

end
