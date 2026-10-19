function [hyp, hasStructure, shapePos, shapeBound, latent] = MF_GP_LearnHyperp(t, y, covFunc, meanFunc, likFunc, infAlg, nfevals, hyp, robust)
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
% hyp,           initial hyperparameter structure (the likelihood and mean parts are used)
% robust,        if true (default: false), fit to convergence (see below)
%
% ---OUTPUTS:
% hyp,           the fitted hyperparameters (NaN if the fit fails)
% hasStructure,  whether the fitted model beats its white-noise special case (always true for the
%                default fit, which does not test); only for robust = true
% shapePos,      positions among the covariance hyperparameters of those that are not identified
%                when hasStructure is false (length scales, periods, and the shape of covRQiso)
% shapeBound,    the values to report for those hyperparameters when hasStructure is false: for a
%                length scale, the sampling interval (a bound: without correlation structure the
%                data do not resolve a length scale above it); NaN for the others
% latent,        the model of the latent function, to make predictions with (fields covFunc, hyp and
%                noisePos): all noise is observation noise, so the covNoise terms are merged with
%                the likelihood noise (see NOTES)
%
% ---NOTES:
% Noise: with a covNoise term in covFunc and a Gaussian likelihood, the model has two noise terms,
% of which only the sum of the variances enters the marginal likelihood (they are fitted equal).
% Predictions at the training inputs would include the covNoise term in the cross-covariance
% (it is a delta function), so that the posterior mean would reproduce part of the noise. For
% predictions, 'latent' therefore drops the covNoise terms from the covariance and merges their
% variances into the likelihood noise: the posterior mean is then that of the smooth component
% (zero for white noise), and the predictive variance is that of a new observation.
%
% Noise floor: the noise standard deviations (the Gaussian likelihood's, and the
% covNoise term's if present) are bounded below by 1% of the standard deviation of
% the data. Without a bound the marginal likelihood of a smooth series is nearly
% flat along a valley in which the noise runs to e^-12 or less, so the fitted noise
% depends on where the optimizer stops and the predictive variance (and ratios with
% it) is set by rounding error. The bounded noise acts as a measurement-noise floor.
% Implemented by clamping the noise hyperparameters inside the objective, with a zero
% gradient for a clamped coordinate.
%
% Robust fit (robust = true): the default fit stops after nfevals function evaluations, which
% is often before the optimizer has converged, so the reported hyperparameters depend on the
% path the optimizer took (and so on round-off) rather than on the data. The robust fit:
% (i) starts from the best of a deterministic grid of length scales (16, log-spaced from the
% sampling interval to half the span) and signal fractions (0.2, 0.6, 0.95) at the total variance of
% the data, scored by the marginal likelihood;
% (ii) minimizes with up to 500 evaluations per restart (up to 3 restarts) of the conjugate-gradient
% minimizer;
% (iii) polishes the optimum with Newton iterations on the gradient (finite-difference Hessian), which
% locate it to round-off (conjugate gradients only locate it to about the square root of machine
% precision over the curvature);
% (iv) bounds the length scales between the sampling interval and the span of the series (a
% shorter or longer length scale cannot be told apart from white noise or a polynomial trend, and
% lets the amplitude and length scale trade off without limit), the shape parameter of
% covRQiso above at alpha = 100 (where it is indistinguishable from covSEiso), and the signal
% amplitudes below at 1% of the standard deviation of the data (as for the noise; a series
% without correlation structure then has a fitted amplitude at this floor, rather than running to zero); and
% (v) tests whether the fitted model beats its white-noise special case (hasStructure): the negative
% log marginal likelihood must be lower than that of zero-mean white noise (variance at its
% maximum-likelihood value) by more than half the 95% chi-squared quantile for the number of
% non-noise kernel hyperparameters. When it does not, the length scales (and shape and period
% parameters) are not identified: a length scale below the sampling interval or a vanishing
% amplitude fits as well, and the likelihood is otherwise flat in them. A length scale is then
% reported as the sampling interval, the bound the data allow (the correlation scale is not
% resolved above one sampling interval); shape and period parameters are NaN. The noise, the
% amplitude (at its floor if the data do not resolve it), and quantities that do not depend on
% the length scale (the marginal likelihood, the predictive mean and variance) are identified,
% and approach those of white noise.
%
% The pre-optimization hyperparameter initialization is data-informed,
% component-by-component, for any covSum of covSEiso/covPeriodic/covRQiso/covNoise;
% unrecognized components fall back to zero.

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
	% Exact Gaussian inference. (Not @infLaplace: with a Gaussian likelihood the
	% Laplace approximation is exact but is computed by Newton iteration, and
	% gpml's infLaplace warm-starts that iteration from a PERSISTENT copy of the
	% previous call's solution, so each fit would depend on whatever series was
	% fitted before it.)
	infAlg = @infGaussLik;
end
if nargin < 7 || isempty(nfevals)
	nfevals = -50; % negative: specifies maximum number of allowed function evaluations
end
% ------------------------------------------------------------------------------

hasStructure = false; latent = [];
% Number of hyperparameters:
s = feval(covFunc{:});
nhps = eval(s);

% Initial values for covariance function:
% Data-informed per-component initialization for any covSum of the component
% covariance functions listed below (an all-zero start can strand the optimizer
% at a substantially worse local optimum).
covFunc1 = covFunc{1};
covFunc2 = covFunc{2};
hyp.cov = zeros(nhps, 1);
noisePos = []; % positions of noise standard deviations among the covariance hyperparameters
alphaPos = []; % positions of the (log) shape parameters of covRQiso
lenPos = []; % positions of the (log) length scales
ampPos = []; % positions of the (log) signal amplitudes
shapePos = []; % positions of the hyperparameters that are not identified without correlation structure
typicalDT = mean(diff(t)); % typical time-step -- used as a length-scale prior
dataSpan = max(t) - min(t);
if strcmp(covFunc1, 'covSum')
	pos = 1;
	for ci = 1:numel(covFunc2)
		if ~ischar(covFunc2{ci})
			% a parameterized component, e.g., {'covMaterniso',3}: gpml reports its number of hyperparameters
			comp = covFunc2{ci};
			numHyp = eval(feval(comp{:}));
			compName = comp{1};
			if isa(compName, 'function_handle')
				compName = func2str(compName);
			end
			if strcmp(compName, 'covMaterniso')
				lenPos = [lenPos, pos]; %#ok<AGROW>
				ampPos = [ampPos, pos + 1]; %#ok<AGROW>
				shapePos = [shapePos, pos]; %#ok<AGROW>
				hyp.cov(pos) = log(typicalDT);   % length-scale
				hyp.cov(pos + 1) = 0;            % log-magnitude
			end
			pos = pos + numHyp;
			continue
		end
		switch covFunc2{ci}
		case 'covSEiso'
			lenPos = [lenPos, pos]; %#ok<AGROW>
			ampPos = [ampPos, pos + 1]; %#ok<AGROW>
			shapePos = [shapePos, pos]; %#ok<AGROW>
			hyp.cov(pos) = log(typicalDT);       % length-scale
			hyp.cov(pos + 1) = 0;                % log-magnitude
			pos = pos + 2;
		case 'covPeriodic'
			shapePos = [shapePos, pos, pos + 1]; %#ok<AGROW>
			ampPos = [ampPos, pos + 2]; %#ok<AGROW>
			hyp.cov(pos) = log(typicalDT);       % length-scale
			hyp.cov(pos + 1) = log(dataSpan / 10); % period (guess: ~10 cycles across the data)
			hyp.cov(pos + 2) = 0;                % log-magnitude
			pos = pos + 3;
		case 'covRQiso'
			lenPos = [lenPos, pos]; %#ok<AGROW>
			ampPos = [ampPos, pos + 1]; %#ok<AGROW>
			hyp.cov(pos) = log(typicalDT);       % length-scale
			hyp.cov(pos + 1) = 0;                % log-magnitude
			hyp.cov(pos + 2) = 0;                % log-alpha (shape)
			alphaPos = [alphaPos, pos + 2]; %#ok<AGROW>
			shapePos = [shapePos, pos, pos + 2]; %#ok<AGROW>
			pos = pos + 3;
		case 'covNoise'
			hyp.cov(pos) = log(0.1);             % noise magnitude
			noisePos = [noisePos, pos]; %#ok<AGROW>
			pos = pos + 1;
		otherwise
			pos = pos + 1; % unrecognized component: leave at zero
		end
	end
end

% ------------------------------------------------------------------------------
% Perform the optimization
% ------------------------------------------------------------------------------
shapeBound = NaN(size(shapePos));
shapeBound(ismember(shapePos, lenPos)) = log(typicalDT); % (the sampling interval, for length scales)
if nargin < 9 || isempty(robust)
	robust = false;
end
if ~robust
	hasStructure = true;
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

	latent = SUB_latent(hyp, covFunc, noisePos);
	return
end

% Robust fit
noiseFloor = log(0.01 * std(y)); % lower bound on log noise standard deviations
bnd = struct('noisePos', noisePos, 'noiseFloor', noiseFloor, 'alphaPos', alphaPos, 'alphaCap', log(100), ...
			'lenPos', lenPos, 'lenLo', log(typicalDT), 'lenHi', log(dataSpan), 'ampPos', ampPos, 'ampFloor', noiseFloor);
hasStructure = false;
hyp.lik = max(hyp.lik, noiseFloor);
hyp.cov(noisePos) = max(hyp.cov(noisePos), noiseFloor);
hyp.cov(alphaPos) = min(hyp.cov(alphaPos), bnd.alphaCap);
hyp.cov(lenPos) = min(max(hyp.cov(lenPos), bnd.lenLo), bnd.lenHi);
hyp.cov(ampPos) = max(hyp.cov(ampPos), bnd.ampFloor);

% Starting point: the best of a grid of length scales and signal fractions, at the total variance of the data
if ~isempty(lenPos) && lenPos(1) == 1 && nhps >= 2
	bestv = Inf; hStart = hyp;
	for logLen = linspace(log(typicalDT), log(dataSpan / 2), 16)
		for phi = [0.2, 0.6, 0.95]
			h = hyp;
			h.cov(1) = logLen;
			h.cov(2) = max(log(sqrt(phi) * std(y)), bnd.ampFloor);
			logNoise = max(log(sqrt((1 - phi) / numel(noisePos)) * std(y)), noiseFloor);
			h.cov(noisePos) = logNoise;
			h.lik = logNoise;
			try
				v = SUB_gpBounded(h, bnd, infAlg, meanFunc, covFunc, likFunc, t, y);
			catch
				v = Inf;
			end
			if isfinite(v) && v < bestv
				bestv = v; hStart = h;
			end
		end
	end
	hyp = hStart;
end

% Minimize, then polish
try
	[hyp, nlZ] = SUB_runOne(hyp, bnd, infAlg, meanFunc, covFunc, likFunc, t, y);
catch emsg
	if strcmp(emsg.identifier, 'MATLAB:nomem')
		error('Not enough memory to fit a Gaussian Process to this data');
	end
	hyp = NaN; return % e.g., a covariance that is not positive definite: not suited to GP fitting
end
if ~isfinite(nlZ)
	hyp = NaN; return
end
try
	hypPolished = SUB_newton(hyp, bnd, infAlg, meanFunc, covFunc, likFunc, t, y);
	nlZ = SUB_gpBounded(hypPolished, bnd, infAlg, meanFunc, covFunc, likFunc, t, y);
	hyp = hypPolished;
catch
	% (keep the unpolished optimum)
end

% Is there any correlation structure to fit? Likelihood-ratio test against the nested white-noise model
n = length(y);
nlZwn = 0.5 * n * log(2 * pi * mean(y.^2)) + 0.5 * n;
numKernel = nhps - numel(noisePos); % number of non-noise kernel hyperparameters
hasStructure = (nlZwn - nlZ) >= gammaincinv(0.95, numKernel / 2); % (half the 95% chi-squared quantile)
latent = SUB_latent(hyp, covFunc, noisePos);

end

% ------------------------------------------------------------------------------
function latent = SUB_latent(hyp, covFunc, noisePos)
	% The model of the latent function: the covNoise terms are observation noise, merged with the
	% likelihood noise (only the sum of their variances matters), so that the covariance of the
	% latent function has no noise term
	latent = struct('covFunc', {covFunc}, 'hyp', hyp, 'noisePos', noisePos);
	if isempty(noisePos) || ~isstruct(hyp) || ~strcmp(covFunc{1}, 'covSum'), return; end
	comps = covFunc{2};
	isNoise = cellfun(@(c) ischar(c) && strcmp(c, 'covNoise'), comps);
	latent.covFunc = {'covSum', comps(~isNoise)};
	latent.hyp.cov = hyp.cov(setdiff(1:numel(hyp.cov), noisePos));
	latent.hyp.lik = 0.5 * log(exp(2 * hyp.lik) + sum(exp(2 * hyp.cov(noisePos))));
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

% ------------------------------------------------------------------------------
function [nlZ, dnlZ] = SUB_gpBounded(hyp, bnd, varargin)
	% gpml's negative log marginal likelihood with the noise hyperparameters clamped
	% at their floor and the covRQiso shape parameters at their cap (and zero gradient
	% wherever they are clamped)
	clampLik = (hyp.lik < bnd.noiseFloor);
	clampCov = bnd.noisePos(hyp.cov(bnd.noisePos) < bnd.noiseFloor);
	atCap = bnd.alphaPos(hyp.cov(bnd.alphaPos) >= bnd.alphaCap);
	hyp.lik = max(hyp.lik, bnd.noiseFloor);
	hyp.cov(bnd.noisePos) = max(hyp.cov(bnd.noisePos), bnd.noiseFloor);
	ampLow = bnd.ampPos(hyp.cov(bnd.ampPos) <= bnd.ampFloor);
	lenLow = bnd.lenPos(hyp.cov(bnd.lenPos) <= bnd.lenLo);
	lenHigh = bnd.lenPos(hyp.cov(bnd.lenPos) >= bnd.lenHi);
	hyp.cov(bnd.alphaPos) = min(hyp.cov(bnd.alphaPos), bnd.alphaCap);
	hyp.cov(bnd.lenPos) = min(max(hyp.cov(bnd.lenPos), bnd.lenLo), bnd.lenHi);
	hyp.cov(bnd.ampPos) = max(hyp.cov(bnd.ampPos), bnd.ampFloor);
	[nlZ, dnlZ] = gp(hyp, varargin{:});
	if clampLik, dnlZ.lik(:) = 0; end
	dnlZ.cov(clampCov) = 0;
	dnlZ.cov(ampLow(dnlZ.cov(ampLow) >= 0)) = 0; % (at the floor, only where the likelihood does not improve above it)
	dnlZ.cov(lenLow(dnlZ.cov(lenLow) >= 0)) = 0;
	dnlZ.cov(lenHigh(dnlZ.cov(lenHigh) <= 0)) = 0;
	dnlZ.cov(atCap(dnlZ.cov(atCap) <= 0)) = 0; % (at the cap, only where the likelihood still improves beyond it)
end

% ------------------------------------------------------------------------------
function [hyp, nlZ] = SUB_runOne(hyp, bnd, infAlg, meanFunc, covFunc, likFunc, t, y)
	% conjugate-gradient minimization of the negative log marginal likelihood, restarted
	% (fresh search direction) until it stops improving, up to 3 times (at most 500 function evaluations each)
	noisePos = bnd.noisePos; noiseFloor = bnd.noiseFloor;
	nlZ = Inf;
	for r = 1:3
		hyp = minimize(hyp, @SUB_gpBounded, -500, bnd, infAlg, meanFunc, covFunc, likFunc, t, y);
		hyp.lik = max(hyp.lik, noiseFloor); % (clamped coordinates can drift: same objective value)
		hyp.cov(noisePos) = max(hyp.cov(noisePos), noiseFloor);
		hyp.cov(bnd.alphaPos) = min(hyp.cov(bnd.alphaPos), bnd.alphaCap);
		hyp.cov(bnd.lenPos) = min(max(hyp.cov(bnd.lenPos), bnd.lenLo), bnd.lenHi);
		hyp.cov(bnd.ampPos) = max(hyp.cov(bnd.ampPos), bnd.ampFloor);
		nlZnew = SUB_gpBounded(hyp, bnd, infAlg, meanFunc, covFunc, likFunc, t, y);
		improved = nlZ - nlZnew;
		nlZ = nlZnew;
		if improved < 1e-10 * max(1, abs(nlZ)), break; end
	end
end

% ------------------------------------------------------------------------------
function hyp = SUB_newton(hyp, bnd, infAlg, meanFunc, covFunc, likFunc, t, y)
	% Newton polish of the conjugate-gradient optimum: conjugate gradients stop when
	% function values can no longer be compared (the hyperparameters are then only
	% known to about the square root of machine precision over the curvature, much worse
	% along a flat valley); Newton iterations on the gradient, which is accurate to round-off,
	% locate the optimum to round-off. The covNoise and Gaussian-likelihood noise terms are
	% interchangeable (only the sum of their variances matters) and stay equal from the symmetric
	% start, so they are tied here; a noise bounded at the floor is held there. Skipped when
	% the Hessian is not positive definite (flat or ridge-like optimum).
	noisePos = bnd.noisePos; noiseFloor = bnd.noiseFloor;
	if numel(noisePos) ~= 1, return; end
	x = hyp.cov(:);
	free = true(size(x));
	if x(noisePos) <= noiseFloor + 1e-9, free(noisePos) = false; end
	free(bnd.alphaPos(x(bnd.alphaPos) >= bnd.alphaCap - 1e-9)) = false; % (shape held at its cap)
	free(bnd.lenPos(x(bnd.lenPos) <= bnd.lenLo + 1e-9 | x(bnd.lenPos) >= bnd.lenHi - 1e-9)) = false;
	free(bnd.ampPos(x(bnd.ampPos) <= bnd.ampFloor + 1e-9)) = false; % (amplitude held at its floor)
	if abs(hyp.lik - x(noisePos)) > 1e-9, return; end % (not tied: leave as is)
	F = find(free);
	for it = 1:20
		[f, g] = SUB_obj(x);
		if max(abs(g(F))) < 1e-10, break; end
		H = zeros(numel(F));
		for a = 1:numel(F)
			h = 1e-5 * max(1, abs(x(F(a))));
			xp = x; xp(F(a)) = xp(F(a)) + h; xm = x; xm(F(a)) = xm(F(a)) - h;
			[~, gp_] = SUB_obj(xp); [~, gm_] = SUB_obj(xm);
			H(:, a) = (gp_(F) - gm_(F)) / (2 * h);
		end
		H = (H + H') / 2;
		ev = eig(H);
		if min(ev) <= 1e-9 * max(ev), break; end % not positive definite: stop
		dx = -H \ g(F);
		m = max(abs(dx));
		if m > 0.5, dx = dx * 0.5 / m; end
		ok = false;
		for ls = 1:10
			xn = x; xn(F) = x(F) + dx;
			fn = SUB_obj(xn);
			if isfinite(fn) && fn <= f + 1e-12 * max(1, abs(f)), ok = true; break; end
			dx = dx / 2;
		end
		if ~ok, break; end
		x = xn;
		if max(abs(dx)) < 1e-12, break; end
	end
	hyp.cov = x; hyp.lik = x(noisePos);

	function [f, g] = SUB_obj(xx)
		hh = hyp; hh.cov = xx; hh.lik = xx(noisePos);
		[f, dd] = SUB_gpBounded(hh, bnd, infAlg, meanFunc, covFunc, likFunc, t, y);
		g = dd.cov(:); g(noisePos) = g(noisePos) + dd.lik;
	end
end
