function out = NL_LargestLyap(y, Nref, maxtstep, past, NNR, embedParams)
% NL_LargestLyap   How quickly initially-close points on the delay-embedded trajectory diverge, as an estimate of the largest Lyapunov exponent.
%
% Rosenstein-style estimate of the largest Lyapunov exponent using TISEAN's
% 'lyap_r' (this operation previously used TSTOOL's 'largelyap'). The series is
% embedded in m dimensions. For each embedded point, the nearest neighbor outside
% a Theiler window (past) is found and the distance between the two is tracked
% as both move forward t = 0,...,maxtstep steps. The divergence curve
% p(t) = <ln d_t> - <ln d_0> (the average log distance, relative to its value at
% t = 0) is, for chaotic dynamics, first linear with slope the largest Lyapunov
% exponent, and then saturates at the size of the attractor. The outputs
% summarize p: its first values, its maximum, the crossings of and times to
% fractions of the maximum, straight-line fits over the rising part (with a
% penalized regression to find the scaling range, attempting to fit as much of
% the range as possible while also achieving the best possible linear fit) and
% an exponential fit.
%
% Both lyap_r and TSTOOL's largelyap track the divergence of each point's single
% nearest neighbor forward in time -- the Rosenstein/Wolf-style construction --
% unlike TISEAN's other Lyapunov tool, 'lyap_k' (Kantz method), which instead
% averages over all neighbors within a range of epsilon-radius neighborhoods.
% Unlike TSTOOL's largelyap, which anchored its output so p(1) = 0, TISEAN's
% lyap_r returns the raw (unanchored) log-divergence; p is re-anchored to
% p(1) = 0 here so that all of the "proportion of max" statistics that follow
% (which assume p rises from ~0) keep working.
%
% ---INPUTS:
% y, the time series to analyze
% Nref, number of randomly-chosen reference points (-1 == all). Accepted for
%       backwards compatibility but no longer affects the computation:
%       'lyap_r' has no equivalent and always uses every valid point as a
%       reference.
% maxtstep, maximum prediction length: {'ac1e', k} for k times the floor of the
%       first 1/e crossing of the autocorrelation function (see BF_GetTau),
%       {'ac', k} for k times its first zero-crossing, or a number of samples
%       (default: {'ac1e',30}). The series must span maxtstep + 2*past <= N/2
%       samples. With {'ac1e', k}, a longer horizon is capped to fit, and the
%       output is NaN only if fewer than max(10, 3 'ac1e' times) remain; with
%       {'ac', k} or a number of samples, the output is NaN instead (capping a
%       zero-crossing horizon lost the correlation of the slopes with known
%       exponents on chaotic flows). {'ac1e', 30} tracked the known exponents of 58
%       dysts flows as well as {'ac', 20} did (Spearman 0.79 vs 0.80 for
%       ve_gradient at N = 1000), and stays defined for finely sampled series,
%       where the zero crossing is long and {'ac', 20} gave NaN for most.
% past, the Theiler window: {'ac1e', k} or {'ac', k} as for maxtstep, or a number
%       of samples (see BF_TheilerWindow; default: {'ac1e',1})
% NNR, number of nearest neighbours. Accepted for backwards compatibility but no
%      longer affects the computation: 'lyap_r' always uses exactly one nearest
%      neighbor per reference point.
% embedParams, input to BF_Embed, how to time-delay-embed the time series, in
%              the form {tau,m}, where string specifiers can indicate standard
%              methods of determining tau ('ac', 'ac1e' or 'mi'; see BF_GetTau) or
%              m ('fnn') (default: {'ac','fnn'})
%
% ---OUTPUTS: statistics of the divergence curve p(t):
% p1, p2, p3, p4, p5: p at steps 0 to 4 (p1 = 0 by construction)
% maxp: the maximum of p
% ncross08max, ncross09max:
%       the number of times p crosses 80% (or 90%) of its maximum
% pcross08max, pcross09max: those numbers as a proportion of the number of steps
% to095max, to09max, to08max, to07max, to05max: the number of steps taken for p
%       to first exceed 95%, 90%, 80%, 70% or 50% of its maximum
% vse_meanabsres, vse_rmsres, vse_gradient, vse_intercept, vse_minbad: a linear
%       fit to p (from 0 to the 95% point of the maximum) with the start and end
%       of the fit varied for the best scaling range: the mean absolute residual,
%       root-mean-square residual, slope, intercept and the minimum value of the
%       penalized error (mean absolute residual minus 0.006 times the number of
%       points)
% ve_meanabsres, ve_rmsres, ve_gradient, ve_intercept, ve_minbad: the same, with
%       only the end of the fit varied (the fit starts at the beginning of p)
% vse_gradient_pertau, ve_gradient_pertau: the two slopes (exponents per sample)
%       multiplied by the (continuous) first 1/e crossing time of the
%       autocorrelation function, i.e. exponents per correlation time. The
%       per-sample slopes scale with the sampling rate; these much less (Lorenz
%       and Rossler flows sampled at 1x to 8x: spread 1.1-2x at N = 5000 and
%       1.2-3.7x at N = 1000, against 4-16x for the per-sample slopes).
% expfit_a, expfit_b, expfit_r2, expfit_adjr2, expfit_rmse: the parameters a and
%       b of a fit p(t) = a*(1 - exp(b*t)) to the whole curve, with its R^2,
%       adjusted R^2 and root-mean-square error
%
% ---REFERENCES:
% A. Wolf et al., "Determining Lyapunov exponents from a time series", Physica D
% 16(3), 285 (1985).
% M. T. Rosenstein, J. J. Collins and C. J. De Luca, "A practical method for
% calculating largest Lyapunov exponents from small data sets", Physica D
% 65(1-2), 117 (1993).
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
%% Check a curve-fitting toolbox license is available:
% ------------------------------------------------------------------------------
BF_CheckToolbox('curve_fitting_toolbox');

doPlot = 0; % whether to plot outputs to a figure

% ------------------------------------------------------------------------------
%% Check inputs, preliminaries:
% ------------------------------------------------------------------------------
N = length(y); % length of time series

% (1) Nref: number of randomly-chosen reference points. No longer used
% (TISEAN's lyap_r always uses every valid point as a reference -- see
% header comment); kept and validated only for call-signature compatibility.
if nargin < 2 || isempty(Nref)
	Nref = 0.5; % use half the length of the time series
end
if Nref < 1 && Nref > 0
	Nref = round(N * Nref); % specify a proportion of time series length
end

% (2) maxtstep: maximum prediction length
if nargin < 3 || isempty(maxtstep)
	maxtstep = {'ac1e', 30}; % (neighbor divergence typically saturates within tens of correlation times)
end
capHorizon = iscell(maxtstep) && strcmp(maxtstep{1}, 'ac1e'); % cap (rather than NaN) a horizon too long for the series
if iscell(maxtstep) % a multiple of the autocorrelation time (as for a Theiler window)
	maxtstep = BF_TheilerWindow(y, maxtstep);
	if isnan(maxtstep)
		warning('No autocorrelation time to set the prediction length')
		out = NaN; return
	end
end
if maxtstep < 1 && maxtstep > 0
	maxtstep = round(N * maxtstep); % specify a proportion of time series length
end
if maxtstep < 10;
	maxtstep = 10; % minimum prediction length; for output stats purposes...
end
if maxtstep > floor(0.5 * N)
	maxtstep = floor(0.5 * N); % can't look further than half the time series length
end

% (3) past/theiler window
if nargin < 4 || isempty(past)
	past = {'ac1e', 1};
end
past = BF_TheilerWindow(y, past, N);
if isnan(past) % no autocorrelation time (e.g. the ACF never decays to 1/e)
	warning('No autocorrelation time to set the Theiler window')
	out = NaN; return
end
if capHorizon && maxtstep + 2 * past > floor(N / 2)
	% Cap an 'ac1e' horizon to fit, keeping at least a few correlation times
	maxtstep = floor(N / 2) - 2 * past;
	if maxtstep < max(10, 3 * BF_GetTau(y, 'ac1e'))
		warning('Time series too short (N = %u) for a prediction horizon with Theiler window %u', N, past);
		out = NaN; return
	end
end
if maxtstep + 2 * past > N / 2
	% Too few correlation times in the series to follow neighbor divergence
	% (lyap_r also becomes very slow as maxtstep + past approaches N). Shortening
	% maxtstep to fit was tested instead: below ~20 autocorrelation times the
	% slope fields lost all correlation with known exponents (dysts flows), so NaN
	warning('Time series too short (N = %u) for maxtstep = %u with Theiler window %u', N, maxtstep, past);
	out = NaN; return
end

% (4) Number of nearest neighbours. No longer used (TISEAN's lyap_r always
% uses exactly one nearest neighbor -- see header comment); kept only for
% call-signature compatibility.
if nargin < 5 || isempty(NNR)
	NNR = 3;
end

% (5) Embedding parameters, embedParams
if nargin < 6 || isempty(embedParams)
	embedParams = {'ac', 'fnn'};
	disp('using default embedding using autocorrelation and cao')
else
	if length(embedParams) ~= 2
		error('Embedding parameters formatted incorrectly -- should be {tau,m}')
	end
end

% ------------------------------------------------------------------------------
%% Resolve the embedding parameters (tau, m)
% ------------------------------------------------------------------------------
tm = BF_Embed(y, embedParams{1}, embedParams{2}, true);
tau = tm(1);
if isnan(tau)
	warning('Could not determine embedding parameters for this time series');
	out = NaN; return
end
m = tm(2);

% ------------------------------------------------------------------------------
%% Run the TISEAN code, lyap_r
% ------------------------------------------------------------------------------
filePath = BF_WriteTempFile(y);
outFilePath = [filePath '.ros'];

[status, res] = BF_TiseanSystem(sprintf('lyap_r -d%u -m%u -t%u -s%u -o %s %s', ...
						  tau, m, past, maxtstep, outFilePath, filePath));

% lyap_r exits 54 when it cannot find a neighbour within its search radius for
% every reference point -- the Rosenstein estimator is not applicable to this
% series (too sparse, or dominated by repeated values / extreme outliers, so
% many embedding vectors are isolated or identical). Treat as a non-result:
if status == 54
	if exist(outFilePath, 'file'), delete(outFilePath); end
	out = NaN; return
end

if isempty(res) || ~isempty(regexp(res, 'command not found', 'once'))
	if exist(outFilePath, 'file'), delete(outFilePath); end
	error('Call to TISEAN function ''lyap_r'' failed.');
end

if ~exist(outFilePath, 'file')
	error('TISEAN function ''lyap_r'' did not produce a .ros output file.');
end

fileInfo = dir(outFilePath);
if fileInfo.bytes == 0
	delete(outFilePath);
	disp('No output obtained from lyap_r');
	out = NaN; return
end

v = dlmread(outFilePath);
delete(outFilePath);
t = v(:, 1);
p = v(:, 2);

if length(p) < 6
	disp('Not enough output from lyap_r to compute statistics');
	out = NaN; return
end

% Re-anchor to p(1) = 0, matching TSTOOL's largelyap convention (see header
% comment) so the "proportion of max" statistics below keep working as before:
p = p - p(1);

if doPlot
	figure('color', 'w');
	box('on');
	plot(t, p, '.-k')
end

% we have the prediction error p as a function of the prediction length...?
% p(tau) = mean over reference points of ln(dist(reference point + tau,
% nearest neighbor + tau) / dist(reference point, nearest neighbor)),
% re-anchored above so p(1) = 0.

% ------------------------------------------------------------------------------
%% Get output stats
% ------------------------------------------------------------------------------

if all(p == 0)
	out = NaN; return
end

% p at lags up to 5
% (note that p1 = 0, so not so useful to record)
for i = 1:5
	% evaluate p(1), p(2), ..., p(5) for the output structure
	out.(sprintf('p%u', i)) = p(i);
end
out.maxp = max(p);

% Number/proportion of crossings at 80% and 90% of maximum
ncrossx = @(x) sum((p(1:end - 1) - x * max(p)) .* (p(2:end) - x * max(p)) < 0);

out.ncross08max = ncrossx(0.8);
out.pcross08max = ncrossx(0.8) / (length(p) - 1);

out.ncross09max = ncrossx(0.9);
out.pcross09max = ncrossx(0.9) / (length(p) - 1);
% out.pcross08max = sum((p(1:end-1)-0.8*max(p)).*(p(2:end)-0.8*max(p)) < 0)/(length(p)-1);
% out.pcross09max = sum((p(1:end-1)-0.9*max(p)).*(p(2:end)-0.9*max(p)) < 0)/(length(p)-1);

% Time taken to get to n% maximum
ttomaxx = @(x) find(p > x * max(p), 1, 'first') - 1;
out.to095max = ttomaxx(0.95);
% out.to095max = find(p > 0.95*max(p),1,'first')-1;
if isempty(out.to095max), out.to095max = NaN; end
out.to09max = ttomaxx(0.9);
% out.to09max = find(p > 0.9*max(p),1,'first')-1;
if isempty(out.to09max), out.to09max = NaN; end
out.to08max = ttomaxx(0.8);
% out.to08max = find(p > 0.8*max(p),1,'first')-1;
if isempty(out.to08max), out.to08max = NaN; end
out.to07max = ttomaxx(0.7);
% out.to07max = find(p > 0.7*max(p),1,'first')-1;
if isempty(out.to07max), out.to07max = NaN; end
out.to05max = ttomaxx(0.5);
% out.to05max = find(p > 0.5*max(p),1,'first')-1;
if isempty(out.to05max), out.to05max = NaN; end

% ------------------------------------------------------------------------------
%% Find scaling region:
% ------------------------------------------------------------------------------
% fit from zero to 95% of maximum...
imax = find(p > 0.95 * max(p), 1, 'first');

if imax <= 3
	% not a suitable range for finding scaling
	% return NaNs for these
	out.vse_meanabsres = NaN;
	out.vse_rmsres = NaN;
	out.vse_gradient = NaN;
	out.vse_intercept = NaN;
	out.vse_minbad = NaN;

	out.ve_meanabsres = NaN;
	out.ve_rmsres = NaN;
	out.ve_gradient = NaN;
	out.ve_intercept = NaN;
	out.ve_minbad = NaN;
else
	t_scal = t(1:imax);
	p_scal = p(1:imax);
	%     pp = polyfit(t_scal,p_scal',1); pfit = pp(1)*t_scal+pp(2);
	% hold on; plot(t_scal,p_scal,'.-r'); hold off
	% hold on; plot(t_scal,pfit,'-r'); hold off;
	% keyboard

	% ------------------------------------------------------------------------------
	%% Adjust start and end times for best scaling
	% ------------------------------------------------------------------------------

	l = imax; % = length(t_scal)
	stptr = 1:floor(l / 2) - 1; % start point must be in the first half (not necessarily, but for here)
	endptr = ceil(l / 2) + 1:l; % end point must be in second half (not necessarily, but for here)
	mybad = zeros(length(stptr), length(endptr));
	for i = 1:length(stptr)
		for j = 1:length(endptr)
			% t_scal/p_scal are both columns here (this used to rely on
			% TSTOOL's spacing() returning t as a row, transposing p_scal to
			% match; with t now a column too, that stray transpose inside
			% lfitbadness turned "pfit - y" into an N-by-N broadcast instead
			% of an N-by-1 residual, so it's dropped):
			mybad(i, j) = lfitbadness(t_scal(stptr(i):endptr(j)), p_scal(stptr(i):endptr(j)));
		end
	end
	[a, b] = find(mybad == min(mybad(:)), 1, 'first'); % this defines the 'best' scaling range (first of any ties)

	% Do the optimum fit again
	t_opt = t_scal(stptr(a):endptr(b));
	p_opt = p_scal(stptr(a):endptr(b));
	pp = polyfit(t_opt, p_opt, 1);
	pfit = pp(1) * t_opt + pp(2);
	res = pfit - p_opt;

	% hold on; plot(t_opt,p_opt,'og'); hold off;
	% hold on; plot(t_opt,pfit,'-g'); hold off;
	% vse == vary start and end times
	out.vse_meanabsres = mean(abs(res));
	out.vse_rmsres = sqrt(mean(res.^2));
	out.vse_gradient = pp(1);
	out.vse_intercept = pp(2);
	out.vse_minbad = min(mybad(:));
	if isempty(out.vse_minbad), out.vse_minbad = NaN; end

	%% Adjust just end time for best scaling
	imin = find(p > 0.50 * max(p), 1, 'first');

	endptr = imin:imax; % end point is at least at 50% mark of maximum
	mybad = zeros(length(endptr), 1);
	for i = 1:length(endptr)
		mybad(i) = lfitbadness(t_scal(1:endptr(i)), p_scal(1:endptr(i)));
	end
	b = find(mybad == min(mybad(:)), 1, 'first'); % this defines the 'best' scaling range (first of any ties)

	% Do the optimum fit again
	t_opt = t_scal(1:endptr(b));
	p_opt = p_scal(1:endptr(b));
	pp = polyfit(t_opt, p_opt, 1);
	pfit = pp(1) * t_opt + pp(2);
	res = pfit - p_opt;

	% hold on; plot(t_opt,p_opt,'om'); hold off;
	% hold on; plot(t_opt,pfit,'-m'); hold off;
	out.ve_meanabsres = mean(abs(res));
	out.ve_rmsres = sqrt(mean(res.^2));
	out.ve_gradient = pp(1);
	out.ve_intercept = pp(2);
	out.ve_minbad = min(mybad(:));
	if isempty(out.ve_minbad), out.ve_minbad = NaN; end

end

% Exponents per correlation time (the continuous first 1/e crossing of the ACF)
tau1e = CO_FirstCrossing(y, 'ac', 1/exp(1), 'continuous');
if isnan(BF_GetTau(y, 'ac1e')), tau1e = NaN; end % the ACF never falls to 1/e
out.vse_gradient_pertau = out.vse_gradient * tau1e;
out.ve_gradient_pertau = out.ve_gradient * tau1e;

% Fit exponential
s = fitoptions('Method', 'NonlinearLeastSquares', 'StartPoint', [max(p) -0.5]);
f = fittype('a*(1-exp(b*x))', 'options', s);
fitWorked = 1;
try
	% t and p are both columns here (this used to rely on TSTOOL's
	% spacing() returning t as a row, transposed back to a column to match
	% fit()'s requirement that X be a column; with t now a column already,
	% that transpose instead turned it into a row, which fit() rejects, so
	% it's dropped):
	[c, gof] = fit(t, p, f);
catch
	fitWorked = 0;
end
if fitWorked
	out.expfit_a = c.a;
	out.expfit_b = c.b;
	out.expfit_r2 = gof.rsquare;
	out.expfit_adjr2 = gof.adjrsquare;
	out.expfit_rmse = gof.rmse;
else
	out.expfit_a = NaN;
	out.expfit_b = NaN;
	out.expfit_r2 = NaN;
	out.expfit_adjr2 = NaN;
	out.expfit_rmse = NaN;
end

if doPlot
	hold('on')
	plot(t, c.a * (1 - exp(c.b * t)), ':r');
	hold('off')
end

% ------------------------------------------------------------------------------
function badness = lfitbadness(x, y, gamma)
	if nargin < 3
		gamma = 0.006; % regularization parameter, gamma, chosen empirically
	end
	pp = polyfit(x, y, 1);
	pfit = pp(1) * x + pp(2);
	res = pfit - y;
	badness = mean(abs(res)) - gamma * length(x); % want to still maximize length(x)
end
% ------------------------------------------------------------------------------

end
