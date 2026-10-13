function out = MF_CompareAR(y, orders, testHow)
% MF_CompareAR   How the out-of-sample error of an AR model changes with its order.
%
% Fits autoregressive (AR) models of a range of orders, and compares the loss of each
% (the mean squared one-step prediction error, see NOTES) when the model fitted to
% a training segment is applied to a test segment. Uses functions from MATLAB's System
% Identification Toolbox: iddata, arxstruc and selstruc. Statistics are taken over the
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
% best_n, the order with the smallest loss (selstruc with criterion 0)
% aic_n, the order that minimizes Akaike's Information Criterion
% bestaic, the minimum value of Akaike's Information Criterion over orders
%
% ---NOTES:
% The loss is the first row of arxstruc's output: the sum of squared one-step
% prediction errors on the test segment, divided by the length of the test segment
% (checked against a least-squares fit by hand, to machine precision). The first
% max(orders) + 1 points of the training and test segments are excluded from the fit
% and from the sum, so that every order is scored on the same points, but the sum is
% still divided by the full test length. The loss is therefore the mean squared
% error scaled by about 1 - (max(orders) + 1)/(test length) (e.g. 0.99 for orders
% 1:10 and 1000 test points), the same for every order.
%
% With testHow = 'all' the models are tested on the data they were trained on, so the
% loss measures in-sample fit: it cannot rise with the model order, and features such
% as minv, propgain1min and where01max mostly describe how fast the fit improves with
% order. Use a training fraction (e.g. 0.5) for a genuine out-of-sample comparison.
%
% If the series is too short for the highest order (the training segment must have more
% than 2*max(orders) + 1 points, and the test segment more than max(orders) + 1), the
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

% ------------------------------------------------------------------------------
%% Check that a System Identification Toolbox license is available:
% ------------------------------------------------------------------------------
BF_CheckToolbox('identification_toolbox');

% Preliminaries
doPlot = false; % set to true to plot outputs
N = length(y); % length of time series, N

%% Check Inputs
% (1) Time series, y
% Convert y to time series object
y = iddata(y, [], 1);

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
% Otherwise arxstruc returns a perfect fit (loss at machine precision, eps) for the high
% orders, or for all of them, and the statistics are an artifact: the output is then NaN.
maxOrder = max(orders(:));
nScoredTrain = size(yTrain, 1) - maxOrder - 1; % points fitted (the first maxOrder + 1 are excluded)
nScoredTest = size(yTest, 1) - maxOrder - 1;
if nScoredTrain <= maxOrder || nScoredTest < 1
	out = NaN; % series too short for this range of orders
	return
end

V = arxstruc(yTrain, yTest, orders);

% ------------------------------------------------------------------------------
%% Output
% ------------------------------------------------------------------------------
% Statistics on V, which contains loss functions at each order (normalized sum of
% squared prediction errors)
v = V(1, 1:end - 1); % the loss function vector over the range of orders

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
%% Use selstruc function to obtain 'best' order measures
% ------------------------------------------------------------------------------
% Get specific 'best' measures
[nn, vmod0] = selstruc(V, 0); % minimizes squared prediction errors
out.best_n = nn;

[nn, vmodaic] = selstruc(V, 'aic'); % minimize Akaike's Information Criterion (AIC)
out.aic_n = nn; % optimum model order minimizing AIC in the range given
% vmodaic is [2 x numOrders]: row 1 holds the AIC of each candidate model, row 2 the orders.
% (Note: the previous form, vmodaic(nn == min(nn)), indexed with a logical scalar and so
%  always returned element 1 -- the AIC of the *first* model, not the best one.)
out.bestaic = min(vmodaic(1, :));

% Using minimum description length is basically the same as using AIC:
% [nn, vmodmdl] = selstruc(V,'mdl'); % minimize Rissanen's Minimum Description Length (MDL)
% out.mdl_n = nn; % optimal model order minimizing MDL in the range given
% out.bestmdl = vmodmdl(nn == min(nn));

end
