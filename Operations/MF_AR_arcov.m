function out = MF_AR_arcov(y,p)
% MF_AR_arcov    Fits an AR model of a given order, p.
%
% Uses arcov code from Matlab's Signal Processing Toolbox.
%
%---INPUTS:
% y, the input time series
% p, the AR model order
%
%---OUTPUTS: include the parameters of the fitted model, the variance estimate
% of a white noise input to the AR model, the root-mean-square (RMS) error of a
% reconstructed time series, and the autocorrelation of residuals.

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

% Does a Signal Processing Toolbox exist?
BF_CheckToolbox('signal_toolbox');

% Check inputs, set defaults:
if nargin < 2 || isempty(p)
    p = 2; % Fit AR(2) model by default
end

%-------------------------------------------------------------------------------
% Fit an AR model using Matlab's Signal Processing Toolbox:
%-------------------------------------------------------------------------------
[a, e] = arcov(y,p);

out.e = e; % variance

% Output fitted parameters up to order, p (+1)
for i = 1:p+1
	out.(sprintf('a%u',i)) = a(i);
end

% ------------------------------------------------------------------------------
%% Residual analysis
% ------------------------------------------------------------------------------
y_est = filter([0, -a(2:end)],1,y);
err = y - y_est; % residuals

out.res_mu = mean(err); % mean error
out.res_std = std(err); % std of error
out.res_AC1 = CO_AutoCorr(err,1,'Fourier'); % autocorrelation of residuals at lag 1
out.res_AC2 = CO_AutoCorr(err,2,'Fourier'); % autocorrelation of residuals at lag 2

end
