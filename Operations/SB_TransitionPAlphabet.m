function out = SB_TransitionPAlphabet(y, numGroups, tau)
% SB_TransitionPAlphabet   How transition probabilities change with alphabet size.
%
% The time series is discretized by quantile separation into numGroups equally
% populated groups, for each alphabet size in a range (default 2 to 10). At each
% size, the matrix T of consecutive-pair frequencies is formed: T(i,j) is the
% number of times group i is followed by group j, divided by N - 1, so T is a
% matrix of joint probabilities that sums to 1 (not a row-normalized transition
% matrix; for equiprobable groups, T(i,j) is about the transition probability
% P(j|i) divided by the alphabet size, as in SB_TransitionMatrix). Six statistics
% of T are computed: the mean and maximum of its diagonal, its trace, its
% asymmetry sum(sum(abs(T - T'))), the trace of its covariance matrix, and the
% standard deviation of its eigenvalues. Exponential decays a*exp(b*x) (and some
% linear fits and change-point statistics) are then fitted to how these change
% with the alphabet size x. Requires the Curve Fitting Toolbox.
%
% ---INPUTS:
% y, the input time series
% numGroups, the number of groups in the coarse-graining: a vector of alphabet
%    sizes to compare across this range (default: 2:10; each must be at least 2).
%    A scalar numGroups, together with a vector tau, is not supported (it errors).
% tau, the time delay at which to analyze the transition matrices (default: 1). We
%    can either downsample the time series at this lag and then do the
%    discretization as normal, or do the discretization and then just look at this
%    discrete lag. Here we do the former (using resample). Can also be a string that
%    sets it from the series: 'ac' (the first zero-crossing of the autocorrelation
%    function, capped at floor(N/50) for a series of length N), 'ac1e' (the floor of
%    its first 1/e crossing), or 'mi' (the smaller of the first minimum of the
%    Kraskov automutual information and the 'ac1e' delay; see BF_GetTau).
%
% ---OUTPUTS:
% A structure with fields (NaN if tau cannot be determined). In the definitions
% below, x is the alphabet size (the entries of numGroups) and T = T_x is the
% joint-probability matrix at that size, so each statistic is a function of x.
% The five fit statistics of each exponential fit a*exp(b*x) are the amplitude a,
% the rate b, R^2, adjusted R^2 and the root-mean-square error of the fit:
% Exponential fit to the mean of the diagonal elements of T, mean_i T(i,i):
% meandiagfexp_a, meandiagfexp_b, meandiagfexp_r2, meandiagfexp_adjr2,
%     meandiagfexp_rmse
% Exponential fit to the maximum of the diagonal elements of T, max_i T(i,i):
% maxdiagfexp_a, maxdiagfexp_b, maxdiagfexp_r2, maxdiagfexp_adjr2,
%     maxdiagfexp_rmse
% Exponential fit to the trace of T, sum_i T(i,i):
% trfexp_a, trfexp_b, trfexp_r2, trfexp_adjr2, trfexp_rmse
% Adjusted R^2 of a linear fit (a*x + b) to the trace of T, over the alphabet sizes
% at which it is above a fifth, or a tenth, of its value for the smallest alphabet
% (NaN if fewer than three sizes qualify):
% trflin5_adjr2, trflin10adjr2
% Asymmetry of T, sum_ij |T(i,j) - T(j,i)|: the slope of a linear fit against
% alphabet size, and the position in the list of alphabet sizes (an index, with 1
% the smallest alphabet size, not the alphabet size itself) of the dividing point
% at which the mean before and after differs most (a t-statistic criterion; NaN if
% the asymmetry is constant or fewer than five sizes are used):
% symd_a, symd_risept
% Trace of the covariance matrix of T, trace(cov(T)): the jump from the first to
% the second alphabet size (value at the second minus value at the first), and the
% exponential fit (excluding the first alphabet size if there is a jump up):
% trcov_jump
% trcovfexp_a, trcovfexp_b, trcovfexp_r2, trcovfexp_adjr2, trcovfexp_rmse
% Exponential fit to the standard deviation of the eigenvalues of T, std(eig(T))
% (the maximum and minimum real eigenvalues of T are no longer fitted):
% stdeigfexp_a, stdeigfexp_b, stdeigfexp_r2, stdeigfexp_adjr2,
%     stdeigfexp_rmse

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
%% Check that a Curve-Fitting Toolbox license is available:
% ------------------------------------------------------------------------------
BF_CheckToolbox('curve_fitting_toolbox');

if nargin < 2 || isempty(numGroups)
	numGroups = (2:10); % compare across alphabet sizes from 2 to 10
end
if nargin < 3 || isempty(tau)
	tau = 1; % use a time-lag of 1
end

N = length(y); % time-series length

if ischar(tau) && ismember(tau, {'ac1e', 'mi'})
	% Adaptive delay: see BF_GetTau
	tau = BF_GetTau(y, tau);
	if isnan(tau)
		out = NaN; return
	end
end
if strcmp(tau, 'ac') % determine tau from first zero of autocorrelation
	tau = CO_FirstCrossing(y, 'ac', 0, 'discrete');
	if isnan(tau)
		out = NaN; return
	end
	if tau > N / 50 % for highly-correlated signals
		tau = floor(N / 50);
	end
end

nfeat = 6; % the number of features calculated at each point
if (length(numGroups) == 1) && (length(tau) > 1) % vary tau
	if numGroups < 2; return; end % need more than 2 groups
	taur = tau; % the tau range
	store = zeros(length(taur), nfeat);

	for i = 1:length(taur)
		tau = taur(i);
		if tau > 1; y = resample(y, 1, tau); end % resample
		yth = SUB_discretize(y, numGroups); % threshold
		store(i, :) = getmeasures(yth);
	end

	error('This setting kind of doesn''t work yet. Sorry.')

elseif (length(tau) == 1) && (length(numGroups) > 1) % vary numGroups
	if min(numGroups) < 2; error('Need more than 2 groups'); end % need more than 2 groups, always
	numGroupsRange = numGroups; % the numGroups range (numGroups is an input vector)
	store = zeros(length(numGroupsRange), nfeat);
	if tau > 1; y = resample(y, 1, tau); end % resample

	for i = 1:length(numGroupsRange)
		numGroups = numGroupsRange(i);
		yth = SUB_discretize(y, numGroups); % thresholded data: yth
		store(i, :) = SUB_getMeasures(yth, numGroups);
	end

	numGroupsRange = numGroupsRange'; % needs to be a column vector for the fitting routines

	% 1) mean of diagonal elements of the transition matrix: shows an exponential
	% decay to zero
	s = fitoptions('Method', 'NonlinearLeastSquares', 'StartPoint', [1, -0.2]);
	f = fittype('a*exp(b*x)', 'options', s);
	[c, gof] = fit(numGroupsRange, store(:, 1), f);
	out.meandiagfexp_a = c.a;
	out.meandiagfexp_b = c.b;
	out.meandiagfexp_r2 = gof.rsquare;
	out.meandiagfexp_adjr2 = gof.adjrsquare;
	out.meandiagfexp_rmse = gof.rmse;

	% 2) maximum of diagonal elements of the transition matrix: shows an exponential
	% decay to zero
	s = fitoptions('Method', 'NonlinearLeastSquares', 'StartPoint', [1, -0.2]);
	f = fittype('a*exp(b*x)', 'options', s);
	[c, gof] = fit(numGroupsRange, store(:, 2), f);
	out.maxdiagfexp_a = c.a;
	out.maxdiagfexp_b = c.b;
	out.maxdiagfexp_r2 = gof.rsquare;
	out.maxdiagfexp_adjr2 = gof.adjrsquare;
	out.maxdiagfexp_rmse = gof.rmse;

	% 3) trace of T
	% fit exponential
	s = fitoptions('Method', 'NonlinearLeastSquares', 'StartPoint', [1, -0.2]);
	f = fittype('a*exp(b*x)', 'options', s);
	[c, gof] = fit(numGroupsRange, store(:, 3), f);
	out.trfexp_a = c.a;
	out.trfexp_b = c.b;
	out.trfexp_r2 = gof.rsquare;
	out.trfexp_adjr2 = gof.adjrsquare;
	out.trfexp_rmse = gof.rmse;

	% Also fit linear from the start to a fifth, a tenth of the starting
	% value
	s = fitoptions('Method', 'NonlinearLeastSquares', 'StartPoint', [-0.05 1]);
	f = fittype('a*x+b', 'options', s);

	r1 = find(store(:, 3) > store(1, 3) / 5);
	if length(r1) > 2
		[~, gof] = fit(numGroupsRange(r1), store(r1, 3), f);
		out.trflin5_adjr2 = gof.adjrsquare;
	else
		out.trflin5_adjr2 = NaN;
	end

	r2 = find(store(:, 3) > store(1, 3) / 10);
	if length(r2) > 2
		[~, gof] = fit(numGroupsRange(r2), store(r2, 3), f);
		out.trflin10adjr2 = gof.adjrsquare;
	else
		out.trflin10adjr2 = NaN;
	end

	% 4) Symmetry; differences in diagonal elements
	% return the slope
	s = fitoptions('Method', 'NonlinearLeastSquares', 'StartPoint', [0.1 0]);
	f = fittype('a*x+b', 'options', s);
	c = fit(numGroupsRange, store(:, 4), f);
	out.symd_a = c.a;

	% return approximately when starts to rise; where means before and
	% after a moving dividing point are most different
	if all(store(:, 4) == store(1, 4)) || length(numGroupsRange) < 5 % all the same, or too few sizes to split
		out.symd_risept = NaN;
	else
		mba = zeros(length(numGroupsRange), 2); % means before and after
		sba = zeros(length(numGroupsRange), 2); % standard deviation before and after
		for i = 3:length(numGroupsRange) - 2
			mba(i, 1) = mean(store(1:i - 1, 4));
			sba(i, 1) = std(store(1:i - 1, 4)) / sqrt(i - 1);
			mba(i, 2) = mean(store(i + 1:end, 4));
			sba(i, 2) = std(store(i + 1:end, 4)) / sqrt(length(numGroupsRange) - i + 1);
		end
		tstats = abs((mba(:, 1) - mba(:, 2)) ./ sqrt(sba(:, 1).^2 + sba(:, 2).^2));
		out.symd_risept = find(tstats == max(tstats), 1, 'first');
	end

	% 5) trace of covariance matrix
	% check jump:
	out.trcov_jump = store(2, 5) - store(1, 5);
	if store(2, 5) > store(1, 5); r1 = 2:length(numGroupsRange); % jump
	else r1 = 1:length(numGroupsRange);
	end
	% fit exponential decay to range without possible first jump
	s = fitoptions('Method', 'NonlinearLeastSquares', 'StartPoint', [1, -0.5]);
	f = fittype('a*exp(b*x)', 'options', s);
	[c, gof] = fit(numGroupsRange(r1), store(r1, 5), f);
	out.trcovfexp_a = c.a;
	out.trcovfexp_b = c.b;
	out.trcovfexp_r2 = gof.rsquare;
	out.trcovfexp_adjr2 = gof.adjrsquare;
	out.trcovfexp_rmse = gof.rmse;

	% 6) Standard deviation of eigenvalues of T
	% Fit an exponential decay
	s = fitoptions('Method', 'NonlinearLeastSquares', 'StartPoint', [1, -0.2]);
	f = fittype('a*exp(b*x)', 'options', s);
	[c, gof] = fit(numGroupsRange, store(:, 6), f);
	out.stdeigfexp_a = c.a;
	out.stdeigfexp_b = c.b;
	out.stdeigfexp_r2 = gof.rsquare;
	out.stdeigfexp_adjr2 = gof.adjrsquare;
	out.stdeigfexp_rmse = gof.rmse;

end

% ------------------------------------------------------------------------------
%% Subfunctions
% ------------------------------------------------------------------------------

function yth = SUB_discretize(y, numGroups)
	% 1) discretize the time series into a number of groups np
	th = quantile(y, linspace(0, 1, numGroups + 1)); % thresholds for dividing the time series values
	th(1) = th(1) - 1; % this ensures the first point is included
	% turn the time series into a set of numbers from 1:numGroups
	yth = zeros(length(y), 1);
	for li = 1:numGroups
		yth(y > th(li) & y <= th(li + 1)) = li;
	end
	if any(yth == 0) % error -- they should all be assigned to a group
		error('Some values were not assigned to a group')
		% yth = []; return;
	end

end

% ------------------------------------------------------------------------------
function out = SUB_getMeasures(yth, numGroups)
	% returns a bunch of metrics on the transition matrix
	N = length(yth);

	% 1) Calculate the one-time transition matrix
	T = zeros(numGroups);
	for j = 1:numGroups
		ri = find(yth == j);
		if isempty(ri) % yth is never j
			T(j, :) = 0;
		else
			if ri(end) == N; ri = ri(1:end - 1); end % looking at next element; remove last point
			for k = 1:numGroups
				T(j, k) = sum(yth(ri + 1) == k); % the next element is of this class
			end
		end
	end
	T = T / (N - 1); % N-1 is appropriate because it's a 1-time transition matrix

	% 2) return some quantities on the transition matrix, T
	%   (i) diagonal elements
	out(1) = mean(diag(T)); % mean of diagonal elements
	out(2) = max(diag(T)); % max of diagonal elements
	out(3) = sum(diag(T)); % sum of diagonal elements (trace)

	%  (ii) measures of symmetry:
	out(4) = sum(sum(abs((T - T')))); % sum of differences of individual elements

	% (iii) measures from covariance matrix:
	out(5) = sum(diag(cov(T))); % trace of covariance matrix

	% (iv) measures from eigenvalues of T
	eigT = eig(T);
	out(6) = std(eigT); % std of eigenvalues

end

end
