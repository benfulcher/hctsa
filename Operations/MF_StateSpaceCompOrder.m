function out = MF_StateSpaceCompOrder(y, maxOrder)
% MF_StateSpaceCompOrder   How the fit of a state-space model improves as its order increases.
%
% Fits state-space models (N4SID subspace identification, BF_SubspaceID) of orders 1,
% 2, ..., maxOrder to the whole time series (all fits are within the sample), and
% returns statistics on how the goodness of fit changes across this range, measured by
% Akaike's information criterion (AIC) and by the loss function (the mean squared
% one-step prediction error). An order at which the model cannot be fitted is left out
% of the summaries (its AIC and loss function are NaN), and the output is NaN only if no
% order can be fitted.
%
% c.f., MF_CompareAR -- does a similar thing for AR models
%
% ---INPUTS:
% y, the input time series
% maxOrder, the maximum model order to consider (default: 10)
%
% ---OUTPUTS:
% minaic: the lowest AIC across orders 1 to maxOrder
% aicopt: the order with the lowest AIC
% minlossfn: the lowest loss function across orders 1 to maxOrder
% lossfnopt: the order with the lowest loss function
% meandiffaic: the mean change in AIC when the order increases by one
% maxdiffaic: the largest increase in AIC when the order increases by one
% mindiffaic: the largest decrease (most negative change) in AIC when the order
%       increases by one
% ndownaic: the number of order increases at which the AIC decreases
% (If some orders cannot be fitted, these are taken over the orders that can; the
% change statistics use only adjacent pairs of orders that both fitted.)
%
% ---NOTES:
% The AIC of order k is log(max(V, 1e-10*mean(y.^2))) + 2*(2k)/N, with V the loss and
% 2k the number of free parameters of an order-k model once the freedom in the choice of
% state coordinates is removed. The floor on V makes orders whose loss is at the level
% of rounding error (an exactly predictable series) tie, so that the lowest of them is
% preferred; lossfnopt is the first order attaining the lowest floored loss.
% meandiffaic is not registered in hctsa's feature set: the mean of the changes in AIC
% telescopes to (AIC(maxOrder) - AIC(1))/(maxOrder - 1) when all orders are fitted, so
% it carries no information beyond the first and last AIC values.

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
if nargin < 2 || isempty(maxOrder)
	maxOrder = 10;
end
y = y(:);
N = length(y);
floorLoss = 1e-10 * mean(y.^2);

% ------------------------------------------------------------------------------
%% Fit the models
% ------------------------------------------------------------------------------
lossfns = NaN(maxOrder, 1); % Loss function
aics = NaN(maxOrder, 1); % Akaike's information criterion
for k = 1:maxOrder
	try
		m = BF_SubspaceID(y, k);
	catch
		continue % this order cannot be fitted
	end
	lossfns(k) = m.loss;
	aics(k) = log(max(m.loss, floorLoss)) + 2 * (2 * k) / N;
end

if all(isnan(aics))
	out = NaN; return % no order could be fitted
end

% ------------------------------------------------------------------------------
%% Outputs
% ------------------------------------------------------------------------------
out.minaic = min(aics);
out.aicopt = find(aics == min(aics), 1, 'first');
out.minlossfn = min(lossfns);
flooredLoss = max(lossfns, floorLoss);
flooredLoss(isnan(lossfns)) = NaN; % (orders that could not be fitted)
out.lossfnopt = find(flooredLoss == min(flooredLoss), 1, 'first');

% Changes in AIC from one order to the next (only adjacent orders that were both fitted)
daics = diff(aics);
daics = daics(~isnan(daics));
if isempty(daics)
	out.meandiffaic = NaN;
	out.maxdiffaic = NaN;
	out.mindiffaic = NaN;
else
	out.meandiffaic = mean(daics);
	out.maxdiffaic = max(daics);
	out.mindiffaic = min(daics);
end
out.ndownaic = sum(daics < 0);

end
