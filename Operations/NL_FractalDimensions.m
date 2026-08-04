function out = NL_FractalDimensions(y, kmin, kmax, Nref, gstart, gend, past, steps, embedParams, randomSeed)
% NL_FractalDimensions    Fractal dimension spectrum, D(q), of a time series.
%
% ---INPUTS:
% y, column vector of time series data
% kmin, minimum number of neighbours for each reference point
% kmax, maximum number of neighbours for each reference point
% Nref, number of randomly-chosen reference points (-1: use all points)
% gstart, starting value for moments
% gend, end value for moments
% past [opt], number of samples to exclude before an after each reference
%             index (default=0)
% steps [opt], number of moments to calculate (default=32);
% embedParams, how to embed the time series using a time-delay reconstruction
% randomSeed [opt], whether (and how) to reset the random seed, using
%             BF_ResetSeed, before choosing reference points (relevant
%             whenever Nref ~= -1, since that involves a random subsample of
%             points). Defaults to 'default' (a fixed seed) so this operation
%             is reproducible by default rather than genuinely stochastic
%             run-to-run.
%
% ---OUTPUTS: include basic statistics of D(q) and q, statistics from a linear fit,
% and an exponential fit of the form D(q) = Aexp(Bq) + C.
%
% Computed natively in MATLAB. This operation previously used TSTOOL's
% 'fracdims', whose actual computation (tstoolbox/@signal/fracdims.m is
% just a thin wrapper; the real work is in the compiled
% mex-dev/GeneralizedDimensionEstimation/gendimest.cpp, vendored in this
% repo under Toolboxes/OpenTSTOOL) turned out to be a well-defined,
% published, citable method rather than an undocumented black box:
% "Generalized Dimensions from Nearest Neighbor Information", P. Schram
% and W. van der Water. Reproduced exactly here:
%
% For each of Nref reference points, find the distances to its 1st..kmax-th
% nearest neighbors (excluding a Theiler window of "past" samples). For
% each moment order gamma swept linearly from gstart to gend (steps
% values), compute the gamma-th moment of the k-th-neighbor distance
% across all reference points, for k = 1:kmax:
%   M(k) = <r_k^gamma>^(1/gamma)   (or exp(<ln r_k>) in the gamma = 0 limit)
% The theoretical relation under an assumed dimension D is
%   M(k) ~ (Gamma(k + gamma/D) / Gamma(k))^(1/gamma)
% (or exp(digamma(k)/D) as gamma -> 0), which needs only its shape -- an
% overall scale factor is fitted separately. D(gamma) is then the value
% that, after re-scaling, best matches this theoretical curve to the
% measured moments M(kmin:kmax) under a robust (log(1+0.5*e^2)) error,
% found via nested 1-D minimization (MATLAB's fminbnd stands in directly
% for the original's hand-rolled Brent's-method minimizer -- the same
% algorithm). Finally q(gamma) = 1 - gamma/D(gamma).
% ------------------------------------------------------------------------------
% Copyright (C) 2020, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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
% Check a curve-fitting toolbox license is available:
% ------------------------------------------------------------------------------
BF_CheckToolbox('curve_fitting_toolbox');

% ------------------------------------------------------------------------------
%% Preliminaries
% ------------------------------------------------------------------------------
N = length(y); % Length of time series
doPlot = 0; % Don't plot results by default

% ------------------------------------------------------------------------------
%% Check inputs
% ------------------------------------------------------------------------------
% (1) Minimum number of neighbours, kmin
if nargin < 2 || isempty(kmin)
	kmin = 3; % default
	fprintf(1, 'Using default, minimum number of neighbours, kmin = %u\n', kmin);
end

% (2) Maximum number of neighbours, kmax
if nargin < 3 || isempty(kmax)
	kmax = 10; % default
	fprintf(1, 'Using default maximum number of neighbours, kmax = %u\n', kmax);
end

% (3) Number of randomly-chosen reference points, Nref
if nargin < 4 || isempty(Nref)
	Nref = 0.2; % default:  20% of the time series length
	fprintf(1, 'Using default number of reference points: Nref = %f\n', Nref);
end
if (Nref > 0) && (Nref < 1)
	Nref = round(N * Nref); % specify a proportion of time series length
end

% (4) moment starting value, gstart
if nargin < 5 || isempty(gstart)
	gstart = 1; % default
	fprintf(1, 'Using default moment starting value, gstart = %u\n', gstart);
end

% (5) moment ending value, gend
if nargin < 6 || isempty(gend)
	gend = 10; % default
	fprintf(1, 'Using default moment ending value, gend = %u\n', gend);
end

% (6) past
if nargin < 7 || isempty(past)
	past = 10; % default
	fprintf(1, 'Using default past correlation exclusion window value, past = %u\n', past);
end

% (7) steps
if nargin < 8 || isempty(steps)
	steps = 32;
end

% (8) Embedding parameters
if nargin < 9 || isempty(embedParams)
	embedParams = {'ac', 'fnn'};
	fprintf(1, 'Using default embedding parameters of autocorrelation for tau and cao method for m\n');
end

% (9) Random seed, for reproducibility of the random reference-point subsample
if nargin < 10 || isempty(randomSeed)
	randomSeed = 'default';
end

% ------------------------------------------------------------------------------
%% Embed the signal (native MATLAB matrix embedding, not a TSTOOL/TISEAN call)
% ------------------------------------------------------------------------------
Y = BF_Embed(y, embedParams{1}, embedParams{2}, false, randomSeed);

if isscalar(Y) && isnan(Y) % embedding failed
	error('Embedding of the %u-sample time series failed', N)
end
[N_embed, m] = size(Y);

if kmax >= N_embed
	% Analogous to TSTOOL's own "too many neighbors ... requested" failure:
	out = NaN; return
end

% ------------------------------------------------------------------------------
%% Resolve the reference points
% ------------------------------------------------------------------------------
if Nref == -1 || Nref >= N_embed
	refIdx = 1:N_embed;
else
	BF_ResetSeed(randomSeed); % for reproducibility of this random subsample
	refIdx = randperm(N_embed, Nref);
end
R = length(refIdx);

% ------------------------------------------------------------------------------
%% For each reference point, find distances to its 1st..kmax-th nearest
%% neighbors (KD-tree, Theiler-window-aware, same approach as
%% NL_localdensity.m/NL_ReturnTime.m)
% ------------------------------------------------------------------------------
kFetch = min(N_embed - 1, kmax + 2 * past + 5);
[idx, dist] = knnsearch(Y, Y(refIdx, :), 'K', kFetch + 1);

distances = zeros(R, kmax); % distances(i,k) = i-th reference point's distance to its k-th nearest neighbor
for ii = 1:R
	i = refIdx(ii);
	validDists = dist(ii, abs(idx(ii, :) - i) > past);
	if length(validDists) < kmax
		allDists = sqrt(sum((Y - Y(i, :)).^2, 2));
		allDists(abs((1:N_embed)' - i) <= past) = Inf;
		validDists = sort(allDists)';
		if sum(isfinite(validDists)) < kmax
			% Not enough valid (Theiler-window-passing) neighbors exist at all:
			out = NaN; return
		end
	end
	distances(ii, :) = validDists(1:kmax);
end

% ------------------------------------------------------------------------------
%% Sweep the moment order and fit a dimension D(gamma) to each
% ------------------------------------------------------------------------------
if (gend - gstart > 0) && (steps > 1)
	gammas = linspace(gstart, gend, steps);
else
	gammas = gstart;
end
numGammas = length(gammas);

Dq = zeros(numGammas, 1);
q = zeros(numGammas, 1);
% Matching gendimest.cpp's own tolerances for the outer (D) and inner
% (scale factor) 1-D minimizations exactly:
outerOpts = optimset('TolX', 1e-4, 'Display', 'off');
innerOpts = optimset('TolX', 1e-5, 'Display', 'off');
for gi = 1:numGammas
	g = gammas(gi);

	% gamma-th moment of the k-th-neighbor distance, across reference points:
	if g == 0
		mom = exp(mean(log(distances), 1));
	else
		mom = mean(distances.^g, 1).^(1 / g);
	end

	DMin = max(0.05, -g / kmin);
	DMax = 128;
	errFcn = @(D) SUB_DimError(D, g, kmin, kmax, mom, innerOpts);
	Dq(gi) = fminbnd(errFcn, DMin, DMax, outerOpts);
	q(gi) = 1 - g / Dq(gi);
end

% ------------------------------------------------------------------------------
% Plot the results in a figure:
% ------------------------------------------------------------------------------
if doPlot
	figure('color', 'w'); box('on');
	plot(q, Dq, 'o-k');
end

% ------------------------------------------------------------------------------
%% Get output statistics:
% ------------------------------------------------------------------------------
out.rangeDq = range(Dq);
out.maxDq = max(Dq);
out.meanDq = mean(Dq);

% out.minq = min(q);
out.maxq = max(q);
out.rangeq = range(q);
out.meanq = mean(q);

% ---Fit linear
% q and Dq are both columns here (this used to rely on TSTOOL's spacing()
% returning q as a row, transposing Dq to match; with q now a column too,
% that stray transpose would turn "p_fit - Dq'" into an N-by-N broadcast
% instead of an N-by-1 residual, so it's dropped):
p = polyfit(q, Dq, 1);
p_fit = q * p(1) + p(2);
res = p_fit - Dq;
out.linfit_a = p(1);
out.linfit_b = p(2);
out.linfit_rmsqres = sqrt(mean(res.^2));

% ---Fit exponential
% s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[range(Dq) -0.5 min(Dq)]);
% f = fittype('a*exp(b*x)+c','options',s);
% [c, gof] = fit(q',Dq,f);
% out.expfit_a = c.a;
% out.expfit_b = c.b;
% out.expfit_c = c.c;
% out.expfit_r2 = gof.rsquare; % I reckon this one is the most important!
% out.expfit_adjr2 = gof.adjrsquare;
% out.expfit_rmse = gof.rmse;

end

% ------------------------------------------------------------------------------
function err = SUB_DimError(D, g, kmin, kmax, mom, fitOpts)
	% Robust curve-fit error between the measured k-th-neighbor-distance
	% moments (mom(kmin:kmax)) and the theoretical curve expected for
	% embedding dimension D, after fitting the theoretical curve's free
	% overall scale factor. Reproduces gendimest.cpp's Error_Function
	% exactly (see this operation's header comment).
	ks = kmin:kmax;
	if g == 0
		z = exp(psi(ks) / D); % psi = digamma function
	else
		z = zeros(size(ks));
		z(1) = 1; % anchored at k = kmin
		running = 1;
		for ii = 1:(length(ks) - 1)
			k = ks(ii);
			running = running * (k + g / D) / k;
			z(ii + 1) = running^(1 / g);
		end
	end

	y = mom(kmin:kmax);

	% Fit the free overall scale factor, a, via the same robust error:
	scaleErrFcn = @(a) sum(log(1 + 0.5 * (y - a * z).^2));
	a = fminbnd(scaleErrFcn, 0, 1e6, fitOpts);

	err = sum(log(1 + 0.5 * (y - a * z).^2));
end
% ------------------------------------------------------------------------------
