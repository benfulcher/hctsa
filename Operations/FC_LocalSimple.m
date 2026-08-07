function out = FC_LocalSimple(y, forecastMeth, trainLength)
% FC_LocalSimple    Simple local time-series forecasting.
%
% Simple predictors using the past trainLength values of the time series to
% predict its next value.
%
% ---INPUTS:
% y, the input time series
%
% forecastMeth, the forecasting method:
%          (i) 'mean': local mean prediction using the past trainLength time-series
%                       values,
%          (ii) 'median': local median prediction using the past trainLength
%                         time-series values
%          (iii) 'lfit': local linear prediction using the past trainLength
%                         time-series values.
%
% trainLength, the number of time-series values to use to forecast the next value
%
% ---OUTPUTS: the mean error, stationarity of residuals, Gaussianity of
% residuals, and their autocorrelation structure.

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
if nargin < 2 || isempty(trainLength)
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
if isempty(evalr)
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
% quantities under the names the rest of the model-fitting family uses. taures itself is
% dropped: it correlates 0.90-0.99 with ac1, whereas its ratio form taurat does not.
residOut = MF_ResidualAnalysis(res, y, 'core');
fields = fieldnames(residOut);
for k = 1:length(fields)
	out.(fields{k}) = residOut.(fields{k});
end

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
