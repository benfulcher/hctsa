function hyp = MF_GP_LearnHyperp(t, y, covFunc, meanFunc, likFunc, infAlg, nfevals, hyp)
% MF_GP_LearnHyperp     Learns Gaussian Process hyperparameters for a time series
%
% Used by main Gaussian Process model fitting operations.
%
% References code 'minimize' from the GAUSSIAN PROCESS REGRESSION AND
% CLASSIFICATION Toolbox version 3.2, which is avilable at:
% http://gaussianprocess.org/gpml/code
%
% ---INPUTS:
%
% t,             time
% y,             data
% covFunc,       the covariance function, formatted as gpml likes it
% meanFunc, the mean function, formatted as gpml likes it
% likFunc, the likelihood function, formatted as gpml likes it
% infAlg, the inference algorithm (in gpml form)
% nfevals,       the number of function evaluations
%
% ---NOTES:
% Noise floor: the noise standard deviations (the Gaussian likelihood's, and the
% covNoise term's if present) are bounded below by 1% of the standard deviation of
% the data. Without a bound the marginal likelihood of a smooth series is nearly
% flat along a valley in which the noise runs to e^-12 or less, so the fitted noise
% depends on where the optimizer stops and the predictive variance (and ratios with
% it) is set by rounding error. The bounded noise acts as a measurement-noise floor.
% Implemented by clamping the noise hyperparameters inside the objective, with a zero
% gradient for a clamped coordinate.
%
% The pre-optimization hyperparameter initialization is data-informed,
% component-by-component, for any covSum of covSEiso/covPeriodic/covRQiso/covNoise;
% unrecognized components fall back to zero. (Data-informed initialization used to be
% done only for the exact covSum {covSEiso,covNoise} pair; every other combination,
% including covSEiso+covPeriodic+covNoise, started from all zeros. This stranded the
% optimizer at a substantially worse local optimum: on real EEG data, refitting the
% periodic combination from a data-informed start changed the negative log marginal
% likelihood by 100s of units, and flipped the periodic-vs-SE-only "which fits better"
% comparison from 2.6% of series to ~33% of series, a far more plausible rate.)

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
	% Exact Gaussian inference. (Was @infLaplace: with a Gaussian likelihood the
	% Laplace approximation is exact but is computed by Newton iteration, and
	% gpml's infLaplace warm-starts that iteration from a PERSISTENT copy of the
	% previous call's solution, so every fit depended on whatever series was
	% fitted before it -- the same series gave outputs differing in the 4th digit
	% from call to call once the optimizer amplified the difference.)
	infAlg = @infGaussLik;
end
if nargin < 7 || isempty(nfevals)
	nfevals = -50; % negative: specifies maximum number of allowed function evaluations
end
% ------------------------------------------------------------------------------

% Number of hyperparameters:
s = feval(covFunc{:});
nhps = eval(s);

% Initial values for covariance function:
% Data-informed per-component initialization for any covSum of the component
% covariance functions listed below. An all-zero start (the previous fallback
% for anything other than the exact covSEiso+covNoise pair) was verified to
% strand the optimizer at a substantially worse local optimum -- e.g. for
% covSEiso+covPeriodic+covNoise, negative log marginal likelihood differed by
% 100s of units between zero- and data-informed starts on real EEG data.
covFunc1 = covFunc{1};
covFunc2 = covFunc{2};
hyp.cov = zeros(nhps, 1);
noisePos = []; % positions of noise standard deviations among the covariance hyperparameters
typicalDT = mean(diff(t)); % typical time-step -- used as a length-scale prior
dataSpan = max(t) - min(t);
if strcmp(covFunc1, 'covSum')
	pos = 1;
	for ci = 1:numel(covFunc2)
		if ~ischar(covFunc2{ci})
			pos = pos + 1; % non-string (e.g., degree-parameterized) component: leave at zero
			continue
		end
		switch covFunc2{ci}
		case 'covSEiso'
			hyp.cov(pos) = log(typicalDT);       % length-scale
			hyp.cov(pos + 1) = 0;                % log-magnitude
			pos = pos + 2;
		case 'covPeriodic'
			hyp.cov(pos) = log(typicalDT);       % length-scale
			hyp.cov(pos + 1) = log(dataSpan / 10); % period (guess: ~10 cycles across the data)
			hyp.cov(pos + 2) = 0;                % log-magnitude
			pos = pos + 3;
		case 'covRQiso'
			hyp.cov(pos) = log(typicalDT);       % length-scale
			hyp.cov(pos + 1) = 0;                % log-magnitude
			hyp.cov(pos + 2) = 0;                % log-alpha (shape)
			pos = pos + 3;
		case 'covNoise'
			hyp.cov(pos) = log(0.1);             % noise magnitude
			noisePos = [noisePos, pos];
			pos = pos + 1;
		otherwise
			pos = pos + 1; % unrecognized component: leave at zero
		end
	end
end

% ------------------------------------------------------------------------------
% Perform the optimization
% ------------------------------------------------------------------------------
try
	% loghyper = minimize(init_loghyper, 'gpr', nfevals, covFunc, t, y);
	noiseFloor = log(0.01 * std(y)); % lower bound on log noise standard deviations
	hyp.lik = max(hyp.lik, noiseFloor);
	hyp.cov(noisePos) = max(hyp.cov(noisePos), noiseFloor);
	hyp = minimize(hyp, @SUB_gpFloored, nfevals, noisePos, noiseFloor, infAlg, meanFunc, covFunc, likFunc, t, y);
	hyp.lik = max(hyp.lik, noiseFloor); % (clamped coordinates can drift: same objective value)
	hyp.cov(noisePos) = max(hyp.cov(noisePos), noiseFloor);
catch emsg
	if strcmp(emsg.identifier, 'MATLAB:posdef')
		fprintf(1, 'Error with lack of positive definite matrix for this function\n');
		hyp = NaN; return % return NaN -- the data is not suited to GP fitting
	elseif strcmp(emsg.identifier, 'MATLAB:nomem')
		error('Not enough memory to fit a Gaussian Process to this data');
	else
		error('Error fitting Gaussian Process to data: %s\n', emsg.message)
	end
end

end

% ------------------------------------------------------------------------------
function [nlZ, dnlZ] = SUB_gpFloored(hyp, noisePos, noiseFloor, varargin)
	% gpml's negative log marginal likelihood with the noise hyperparameters clamped
	% at noiseFloor (and zero gradient wherever they are clamped)
	clampLik = (hyp.lik < noiseFloor);
	clampCov = noisePos(hyp.cov(noisePos) < noiseFloor);
	hyp.lik = max(hyp.lik, noiseFloor);
	hyp.cov(noisePos) = max(hyp.cov(noisePos), noiseFloor);
	[nlZ, dnlZ] = gp(hyp, varargin{:});
	if clampLik, dnlZ.lik(:) = 0; end
	dnlZ.cov(clampCov) = 0;
end
