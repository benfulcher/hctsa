% MF_AR_arcov
% 
% Fits an AR model of a given order, p, using arcov code from Matlab's Signal
% Processing Toolbox.
% 
% INPUTS:
% y, the input time series
% p, the AR model order
% 
% Outputs include the parameters of the fitted model, the variance estimate of a
% white noise input to the AR model, the root-mean-square (RMS) error of a
% reconstructed time series, and the autocorrelation of residuals.
% 
% ------------------------------------------------------------------------------
% Copyright (C) 2013,  Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones., "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2010). DOI: 10.1098/rsif.2013.0048
%
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
% details.
% 
% You should have received a copy of the GNU General Public License along with
% this program.  If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

function out = MF_AR_arcov(y,p)
% Ben Fulcher, 2009

N = length(y); % Length of input time series

if nargin < 2 || isempty(p)
    p = 2; % Fit AR(2) model by default
end

% Fit an AR model using Matlab's Signal Processing Toolbox:
[a, e] = arcov(y,p);

out.e = e; % variance

% Output fitted parameters up to order, p (+1)
for i = 1:p+1
	eval(sprintf('out.a%u = a(%u);',i,i));
end

%% Residual analysis
y_est = filter([0, -a(2:end)],1,y);
err = y - y_est; % residuals

out.rms = sqrt(sum(err.^2)/N); % RMS error
out.mu = mean(err); % mean error
out.std = std(err); % std of error
out.AC1 = CO_AutoCorr(err,1); % autocorrelation of residuals at lag 1
out.AC2 = CO_AutoCorr(err,2); % autocorrelation of residuals at lag 2

end