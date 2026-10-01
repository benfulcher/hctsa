function out = FC_LocalSimple(y, forecastMeth, trainLength)
% FC_LocalSimple   How well a simple local rule forecasts the next value.
%
% Predicts each value of the time series from the trainLength values just before
% it, using a simple rule (the mean, the median, or a straight line fitted to those
% values and extended one step). The residuals are the prediction minus the data
% (e = yp - y), and the outputs describe them with the shared residual summary
% MF_ResidualAnalysis (at its 'core' level), plus the Gaussianity of their
% distribution. The first trainLength values are used for training only and are not
% forecast. If the series is too short to forecast, a NaN is returned.
%
% ---INPUTS:
% y, the input time series
%
% forecastMeth, the forecasting method:
%          (i) 'mean': the mean of the past trainLength values (default),
%          (ii) 'median': the median of the past trainLength values,
%          (iii) 'lfit': the next point on a straight line fitted to the past
%                         trainLength values.
%
% trainLength, the number of past values used to forecast the next value (default
%          3), or 'ac' to use the first zero-crossing of the autocorrelation
%          function of y (discrete, from CO_FirstCrossing).
%
% ---OUTPUTS:
% meane, mean of the residuals (the bias of the forecast)
% meanabs, mean absolute residual
% stde, standard deviation of the residuals
% maxonstd, largest absolute residual, in units of the residual standard deviation
% ac1, ac2, ac3: autocorrelation of the (z-scored) residuals at lags 1, 2 and 3
% propbth, proportion of the residual autocorrelations at lags 1 to 25 within the
%          significance band +/- 2.6/sqrt(N)
% taurat, decorrelation time of the residuals (first zero-crossing of their
%          autocorrelation function) divided by that of the time series
% sws, standard deviation across 5 windows of the local standard deviation of the
%          residuals, relative to their overall standard deviation
% swm, standard deviation across 5 windows of the local mean of the residuals,
%          relative to their overall standard deviation
% normr2, R^2 of a Gaussian fit (DN_SimpleFit) to the kernel-smoothed distribution
%          of the residuals: a measure of their Gaussianity

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
% Check inputs
% ------------------------------------------------------------------------------
% Forecasting method, forecastMeth
if nargin < 2 || isempty(forecastMeth)
	forecastMeth = 'mean';
end
% Number of samples to train with, trainLength
if nargin < 3 || isempty(trainLength)
	trainLength = 3;
end

N = length(y); % Time-series length

% ------------------------------------------------------------------------------
% Do the local prediction
% ------------------------------------------------------------------------------
if strcmp(trainLength, 'ac')
	% Make it first zero-crossing of ACF:
	lp = CO_FirstCrossing(y, 'ac', 0, 'discrete');
else
	lp = trainLength; % the length of the subsegment preceeding to use to predict the subsequent value
end
evalr = lp + 1:N; % range over which to evaluate the forecast
if isnan(lp) || isempty(evalr)
	warning('This time series is too short for forecasting');
	out = NaN;
	return
end
res = zeros(length(evalr), 1); % residuals

switch forecastMeth
	case 'mean'
		for i = 1:length(evalr)
			res(i) = mean(y(evalr(i) - lp:evalr(i) - 1)) - y(evalr(i)); % prediction - value
		end
	case 'median'
		for i = 1:length(evalr)
			res(i) = median(y(evalr(i) - lp:evalr(i) - 1)) - y(evalr(i)); % prediction - value
		end
	case 'lfit'
		for i = 1:length(evalr)
			% Fit linear
			warning('off', 'MATLAB:polyfit:PolyNotUnique'); % Disable (potentially important ;)) warning
			p = polyfit((1:lp)', y(evalr(i) - lp:evalr(i) - 1), 1);
			warning('on', 'MATLAB:polyfit:PolyNotUnique'); % Re-enable warning
			res(i) = polyval(p, lp + 1) - y(evalr(i)); % prediction - value
		end
	otherwise
		error('Unknown forecasting method ''%s''', forecastMeth);
end

% out=res;
% plot(res);

% ------------------------------------------------------------------------------
% Output statistics on the residuals, res
% ------------------------------------------------------------------------------

% Report the residuals through the shared contract, at the cheap 'core' level. This
% replaces the hand-rolled meanerr, stderr, meanabserr, sws, swm, ac1, ac2, taures and
% tauresrat, which are now meane, stde, meanabs, sws, swm, ac1, ac2 and taurat -- the same
% quantities under the names the rest of the model-fitting family uses.
residOut = MF_ResidualAnalysis(res, y, 'core');
fields = fieldnames(residOut);
for k = 1:length(fields)
	out.(fields{k}) = residOut.(fields{k});
end

% (Dropped: taures, the residual decorrelation time -- it correlates 0.90-0.99 with ac1.
%  Its ratio form taurat, which does not, is retained via the shared contract above.)

% Normality of the residuals, as the r-squared of a Gaussian fit to their distribution.
% (Renamed from gofr2, which read as a goodness-of-fit measure for the *forecast*; it is
% not -- it is a distributional statistic, so it sits alongside the contract's normksstat
% rather than alongside stde.)
tmp = DN_SimpleFit(res, 'gauss1', 0);
if ~isstruct(tmp) && isnan(tmp) % fitting failed
	out.normr2 = NaN;
else
	out.normr2 = tmp.r2; % r-squared
end

end
