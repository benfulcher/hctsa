function out = MF_AR_arcov(y, p)
% MF_AR_arcov   An AR model of a given order fitted to the time series.
%
% Fits an autoregressive (AR) model of order p, x(t) + a(2)*x(t-1) + ... +
% a(p+1)*x(t-p) = e(t), by least-squares fitting of the one-step prediction (the
% covariance method, arcov from MATLAB's Signal Processing Toolbox). The outputs are
% the fitted polynomial coefficients, the variance of the white noise that drives
% the model, and statistics of the residuals (the data minus the one-step
% prediction), from the shared residual summary MF_ResidualAnalysis ('core' level).
%
% ---INPUTS:
% y, the input time series
% p, the AR model order (default 2)
%
% ---OUTPUTS:
% noisevar, the variance of the white noise input to the fitted AR model
% a2, a3, a4, a5, a6 (up to a(p+1)): the fitted AR polynomial coefficients; a(k+1)
%       is the negative of the usual AR coefficient on the lag-k value. (a1 is
%       always 1 and is also returned.)
% meane, mean of the residuals (note the sign: data minus prediction)
% meanabs, mean absolute residual
% stde, standard deviation of the residuals
% maxonstd, largest absolute residual, in units of the residual standard deviation
% ac1, ac2, ac3: autocorrelation of the (z-scored) residuals at lags 1, 2 and 3
% propbth, proportion of the residual autocorrelations at lags 1 to 25 within the
%       significance band +/- 2.6/sqrt(N)
% taurat, decorrelation time of the residuals (first zero-crossing of their
%       autocorrelation function) divided by that of the time series
% sws, standard deviation across 5 windows of the local standard deviation of the
%       residuals, relative to their overall standard deviation
% swm, standard deviation across 5 windows of the local mean of the residuals,
%       relative to their overall standard deviation
%
% ---NOTES:
% The first p residuals are computed with the unseen values before the start of the
% series set to zero, and are not excluded from the residual statistics.

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

% Does a Signal Processing Toolbox exist?
BF_CheckToolbox('signal_toolbox');

% Check inputs, set defaults:
if nargin < 2 || isempty(p)
	p = 2; % Fit AR(2) model by default
end

% -------------------------------------------------------------------------------
% Fit an AR model using Matlab's Signal Processing Toolbox:
% -------------------------------------------------------------------------------
[a, e] = arcov(y, p);

% Variance of the white noise driving the fitted AR process. Named noisevar to match
% MF_armax and MF_StateSpace_n4sid, which report the same quantity:
out.noisevar = e;

% Output fitted parameters up to order, p (+1)
for i = 1:p + 1
	out.(sprintf('a%u', i)) = a(i);
end

% ------------------------------------------------------------------------------
%% Residual analysis
% ------------------------------------------------------------------------------
y_est = filter([0, -a(2:end)], 1, y);
err = y - y_est; % residuals

% Report the residuals through the shared contract, at the cheap 'core' level (this
% operation is meant to stay cheap). Replaces the four hand-rolled statistics res_mu,
% res_std, res_AC1 and res_AC2, which are now meane, stde, ac1 and ac2.
residOut = MF_ResidualAnalysis(err, y, 'core');
fields = fieldnames(residOut);
for k = 1:length(fields)
	out.(fields{k}) = residOut.(fields{k});
end

end
