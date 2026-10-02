function out = MF_StateSpaceCompOrder(y, maxOrder)
% MF_StateSpaceCompOrder   How the fit of a state-space model improves as its order increases.
%
% Fits state space models using n4sid (from Matlab's System Identification
% Toolbox) of orders 1, 2, ..., maxOrder to the whole time series (all fits are
% within the sample), and returns statistics on how the goodness of fit changes
% across this range, measured by Akaike's information criterion (AIC) and by the
% loss function (the estimated variance of the one-step prediction error).
% An order at which the model cannot be fitted is left out of the summaries (its
% AIC and loss function are NaN), and the output is NaN only if no order can be fitted.
%
% c.f., MF_CompareAR -- does a similar thing for AR models
% Uses the functions iddata, n4sid, and aic from Matlab's System Identification
% Toolbox
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
% Akaike's final prediction error is also computed at each order but is not output.

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

% ------------------------------------------------------------------------------
% Check inputs:
% ------------------------------------------------------------------------------
% Maximum model order, maxOrder (compare models from order 1 up to
%           this)
if nargin < 2 || isempty(maxOrder)
	maxOrder = 10;
end
% orders = 1:maxOrder;

% ------------------------------------------------------------------------------
%% Preliminaries
% ------------------------------------------------------------------------------

% Convert y to time series object
y = iddata(y, [], 1);

% ------------------------------------------------------------------------------
%% Fit the state space models, returning basic fit statistics as we go
% ------------------------------------------------------------------------------

% Initialize statistics -- all within-sample statistics. Could also fit on
% a portion and then predict on another...

% noisevars = zeros(maxOrder,1); % Noise variance -- for us the same as
% loss fn
lossfns = NaN(maxOrder, 1); % Loss function
fpes = NaN(maxOrder, 1); % Akaike's final prediction error
aics = NaN(maxOrder, 1); % Akaike's information criterion

for k = 1:maxOrder
	% Fit the state space model for this order, k
	try
		m = n4sid(y, k);
	catch
		% Data-dependent (n4sid could not fit this series at this order), so NaN
		% at this order only, rather than error(), per the NaN-vs-error convention:
		warning('State-space model fitting failed for k = %u', k);
		continue
	end

	lossfns(k) = m.EstimationInfo.LossFcn;
	fpes(k) = m.EstimationInfo.FPE;
	aics(k) = aic(m);
end

if all(isnan(aics))
	out = NaN; return % no order could be fitted
end

% Optimum model orders (over the orders that could be fitted)
out.minaic = min(aics);
out.aicopt = find(aics == min(aics), 1, 'first');
% out.minbic = min(bics);
% out.bicopt = find(bics == min(bics), 1, 'first');
out.minlossfn = min(lossfns);
out.lossfnopt = find(lossfns == min(lossfns), 1, 'first');

% Curve change summary statistics (over pairs of adjacent orders that both fitted)
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
