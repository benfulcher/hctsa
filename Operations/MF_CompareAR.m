function out = MF_CompareAR(y, orders, testHow)
% MF_CompareAR   How the out-of-sample error of an AR model changes with its order.
%
% Fits autoregressive (AR) models of a range of orders, and compares the loss of each
% (the mean squared one-step prediction error, see NOTES) when the model fitted to
% a training segment is applied to a test segment. Statistics are taken over the
% loss as a function of model order, v.
%
% ---INPUTS:
% y, vector of time-series data
%
% orders, a vector of model orders to compare (default 1:10)
%
% testHow, a fraction of the time series to train on (the model is tested on the
%          remaining portion), or the string 'all' to train and test on all the
%          data (default). With 'all' the loss measures in-sample fit, not
%          out-of-sample prediction.
%
% ---OUTPUTS:
% maxv, minv, meanv, medianv: the maximum, minimum, mean and median of the loss over
%       orders
% propgain1min, the proportion of the first order's loss removed by the best order,
%       1 - min(v)/v(1) (between 0 and 1; 1 for a perfectly predictable series; NaN if
%       the first loss is zero or not finite)
% medonmax, the median loss divided by the maximum loss (in (0, 1], approaching 1 when
%       the loss is insensitive to the order; NaN if the maximum is zero or not finite)
% meandiff, stddiff, maxdiff, meddiff: the mean, standard deviation, maximum absolute
%       value and median of the change in loss from one order to the next
% minstdfromi, the minimum (over starting orders i) of the standard error of the loss
%       over orders i onward, std(v(i:end))/sqrt(length(v)-i+1), ignoring zeros
% where01max, the first position in the list of orders from which that standard error
%       is below 10% of its maximum (NaN if none)
% whereen4, the first position in the list of orders from which it is below 1e-4
%       (NaN if none)
% best_n, the order with the smallest loss
% aic_n, the order that minimizes Akaike's Information Criterion,
%       log(loss*(1 + 2*order/n)) with n the number of test points scored
% bestaic, the minimum value of that criterion over orders
%
% ---NOTES:
% For each order p, an AR(p) model without a constant is fitted by least squares to the
% training segment and scored on the test segment. All orders are fitted to and scored
% on the same points: the first max(orders) points of each segment are excluded (they
% lack a full set of max(orders) past values). The loss is the mean squared one-step
% prediction error over the n scored test points.
%
% With testHow = 'all' the models are tested on the data they were trained on, so the
% loss measures in-sample fit: it cannot rise with the model order, and features such
% as minv, propgain1min and where01max mostly describe how fast the fit improves with
% order. Use a training fraction (e.g. 0.5) for a genuine out-of-sample comparison.
%
% best_n and aic_n (AIC(p) = log(max(v(p), floor)*(1 + 2p/n))) are chosen from the
% losses floored at 1e-10 times the mean square of the test segment, so that orders
% whose loss is at the level of rounding error (an exactly predictable series) tie and
% the lowest of them is chosen; bestaic is the minimum of that criterion.
%
% If the series is too short for the highest order (the training segment must have more
% than 2*max(orders) points, and the test segment more than max(orders)), the
% highest-order models interpolate the training data and the loss is at machine precision;
% NaN is returned for every output.

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

% Preliminaries
doPlot = false; % set to true to plot outputs
N = length(y); % length of time series, N

%% Check Inputs
% (1) Time series, y
y = y(:);

% (2) Model orders:
if nargin < 2 || isempty(orders)
	orders = (1:10)';
end
if size(orders, 1) == 1
	orders = orders'; % make sure a column vector
end

% (3) testHow -- either all (trains and tests on the whole time series);
% or a proportion of the time series to train on; will test on the
% remaining portion.
if nargin < 3 || isempty(testHow)
	testHow = 'all';
end

% ------------------------------------------------------------------------------
%% Run
% ------------------------------------------------------------------------------
% Get normalized prediction errors, V, from training to test set for each
% model order
% This could be done for model residuals using code in, say,
% MF_StateSpaceCompOrder, or MF_StateSpace_n4sid...

if ischar(testHow)
	if strcmp(testHow, 'all')
		yTrain = y;
		yTest = y;
	else
		error('Unknown testing set specifier ''%s''', testHow);
	end
else
	% use first <proportion> to train, rest to test
	co = floor(N * testHow); % cutoff
	yTrain = y(1:co);
	yTest = y(co + 1:end);
end

% The loss is only meaningful if the highest-order model is identifiable from the training
% segment (more points fitted than parameters) and the test segment has points to score.
% Otherwise the highest orders interpolate the training data (loss at machine precision)
% and the statistics are an artifact: the output is then NaN.
maxOrder = max(orders(:));
nScoredTrain = length(yTrain) - maxOrder; % points fitted (the first maxOrder are excluded)
nScoredTest = length(yTest) - maxOrder;
if nScoredTrain <= maxOrder || nScoredTest < 1
	out = NaN; % series too short for this range of orders
	return
end

% Least-squares AR fits (targets from maxOrder + 1) scored on the test segment
numOrders = length(orders);
v = zeros(1, numOrders);
for k = 1:numOrders
	p = orders(k);
	tTrain = (maxOrder + 1:length(yTrain))';
	Xtrain = zeros(length(tTrain), p);
	tTest = (maxOrder + 1:length(yTest))';
	Xtest = zeros(length(tTest), p);
	for j = 1:p
		Xtrain(:, j) = yTrain(tTrain - j);
		Xtest(:, j) = yTest(tTest - j);
	end
	phi = BF_Pinv(Xtrain) * yTrain(tTrain);
	err = yTest(tTest) - Xtest * phi;
	v(k) = mean(err.^2);
end
nTest = nScoredTest;

% ------------------------------------------------------------------------------
%% Output
% ------------------------------------------------------------------------------
% Statistics on V, which contains loss functions at each order (normalized sum of
% squared prediction errors)

out.maxv = max(v);
out.minv = min(v);
out.meanv = mean(v);
out.medianv = median(v);
% Bounded forms of the two loss ratios: the proportional improvement from the first to the
% best order, and the median over the maximum (a zero or non-finite denominator is NaN)
if v(1) > 0 && isfinite(v(1))
	out.propgain1min = 1 - min(v) / v(1);
else
	out.propgain1min = NaN;
end
if max(v) > 0 && isfinite(max(v))
	out.medonmax = median(v) / max(v);
else
	out.medonmax = NaN;
end
out.meandiff = mean(diff(v));
out.stddiff = std(diff(v));
out.maxdiff = max(abs(diff(v)));
out.meddiff = median(diff(v));

% where does it steady off?
stdfromi = zeros(length(v), 1);
for i = 1:length(stdfromi)
	stdfromi(i) = std(v(i:end)) / sqrt(length(v) - i + 1);
end
out.minstdfromi = min(stdfromi(stdfromi > 0));
if isempty(out.minstdfromi), out.minstdfromi = NaN; end
out.where01max = find(stdfromi < max(stdfromi) * 0.1, 1, 'first');
if isempty(out.where01max), out.where01max = NaN; end
out.whereen4 = find(stdfromi < 1e-4, 1, 'first');
if isempty(out.whereen4), out.whereen4 = NaN; end

% ------------------------------------------------------------------------------
%% Plotting
% ------------------------------------------------------------------------------
if doPlot
	plot(v);
	plot(stdfromi, 'r');
end

% ------------------------------------------------------------------------------
%% 'Best' orders by loss and by AIC
% ------------------------------------------------------------------------------
% Losses floored at the level of rounding error, so that exactly predictable series tie
floorV = 1e-10 * mean(yTest.^2);
vFloored = max(v, floorV);
[~, ix] = min(vFloored);
out.best_n = orders(ix);
aic = log(vFloored .* (1 + 2 * orders' / nTest));
[out.bestaic, ix] = min(aic);
out.aic_n = orders(ix);

end
