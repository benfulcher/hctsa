function out = MF_StateSpaceCompOrder(y,maxOrder)
% MF_StateSpaceCompOrder    Change in goodness of fit across different state space models.
%
% Fits state space models using n4sid (from Matlab's System Identification
% Toolbox) of orders 1, 2, ..., maxOrder and returns statistics on how the
% goodness of fit changes across this range.
%
% c.f., MF_CompareAR -- does a similar thing for AR models
% Uses the functions iddata, n4sid, and aic from Matlab's System Identification
% Toolbox
%
%---INPUTS:
% y, the input time series
% maxOrder, the maximum model order to consider.

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

% ------------------------------------------------------------------------------
%% Check that a System Identification Toolbox license is available:
% ------------------------------------------------------------------------------
BF_CheckToolbox('identification_toolbox')

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
y = iddata(y,[],1);

% ------------------------------------------------------------------------------
%% Fit the state space models, returning basic fit statistics as we go
% ------------------------------------------------------------------------------

% Initialize statistics -- all within-sample statistics. Could also fit on
% a portion and then predict on another...

% noisevars = zeros(maxOrder,1); % Noise variance -- for us the same as
% loss fn
lossfns = zeros(maxOrder,1); % Loss function
fpes = zeros(maxOrder,1); % Akaike's final prediction error
aics = zeros(maxOrder,1); % Akaike's information criterion

for k = 1:maxOrder
    % Fit the state space model for this order, k
    try
        m = n4sid(y,k);
    catch emsg
        error('Model fitting failed for k = %u',k)
    end

    lossfns(k) = m.EstimationInfo.LossFcn;
    fpes(k) = m.EstimationInfo.FPE;
    aics(k) = aic(m);
end

% Optimum model orders
out.minaic = min(aics);
out.aicopt = find(aics == min(aics), 1, 'first');
% out.minbic = min(bics);
% out.bicopt = find(bics == min(bics), 1, 'first');
out.minlossfn = min(lossfns);
out.lossfnopt = find(lossfns == min(lossfns), 1, 'first');

% Curve change summary statistics
out.meandiffaic = mean(diff(aics));
out.maxdiffaic = max(diff(aics));
out.mindiffaic = min(diff(aics));
out.ndownaic = sum(diff(aics) < 0);

end
