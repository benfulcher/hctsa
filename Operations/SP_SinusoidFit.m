function out = SP_SinusoidFit(y, model)
% SP_SinusoidFit   Fit sinusoids or a Fourier series to the time series.
%
% Fits a sum of 1-3 sinusoids, or a Fourier series with 1-3 terms, to the time
% series y as a function of its time index t = 1:N (the values are fitted in the
% order in which they occur, so the result depends on the temporal ordering of
% the data). Uses the 'fit' function from Matlab's Curve Fitting Toolbox.
%
% The fitted models are:
%   'sinK':     a sum of K sinusoids, sum_{i=1..K} a_i*sin(2*pi*f_i*t + c_i), with
%               free amplitudes a_i, phases c_i, and frequencies f_i in
%               [1/N, 1/2 - 1/(2N)] cycles per sample,
%   'fourierK': a K-term Fourier series,
%               a0 + sum_{i=1..K} (a_i*cos(i*w*t) + b_i*sin(i*w*t)),
%               with a single fitted fundamental frequency w (so that the terms
%               are harmonically related).
%
% Goodness of fit is summarized using the root mean square error (and R^2), and
% the residuals of the fit are then characterized by their autocorrelation at
% lags 1 and 2 and by a runs test, which together reveal remaining temporal
% structure that the model has not captured (e.g., whether a periodic component
% has been fully explained).
%
% ---INPUTS:
% y, the input time series (a vector; row vectors are converted to columns)
% model, the model to fit:
%       (i) 'sin1': a single sinusoid
%       (ii) 'sin2': a sum of two sinusoids
%       (iii) 'sin3': a sum of three sinusoids
%       (iv) 'fourier1': a Fourier series with one term
%       (v) 'fourier2': a Fourier series with two terms
%       (vi) 'fourier3': a Fourier series with three terms
%
% ---OUTPUTS: a structure containing
% r2, the R^2 of the fit, relative to the mean of the series. The sinK models have
%         no intercept (constant) term, so r2 can be negative for them
% adjr2, the degrees-of-freedom-adjusted R^2
% rmse, the root mean square error of the fit (the residual sum of squares divided
%         by the degrees of freedom of the error, as in the Curve Fitting Toolbox)
% resAC1, the autocorrelation of the residuals at lag 1 (using the 'Fourier'
%         method of CO_AutoCorr)
% resAC2, the autocorrelation of the residuals at lag 2
% resrunsz, the signed z-statistic of a runs test on the residuals (BF_RunsZ):
%         negative when the residuals have fewer runs about their median than
%         expected for a random order (slowly varying residuals)
% resAC1, resAC2 and resrunsz are NaN if the residuals have no variation beyond
% floating-point round-off, as for an exact fit (BF_ResidualStats).
% If the model cannot be fitted (NaN or Inf computed by the model function, or
% fewer than 3K+1 samples for K sinusoids), NaN is returned instead of a structure.
%
% ---NOTES:
% These time-series models are kept separate from DN_SimpleFit, which fits the
% distribution of values: that distribution is unaffected by temporal ordering,
% whereas these fits are not. r2 and adjr2 are not registered for the
% sin1/sin2/sin3 mops (rmse is).
% The sinusoid frequencies are bounded below by 1/N cycles per sample, so that each
% sinusoid completes at least one full cycle within the series: these models measure
% how well the series is described by oscillations, and a component with less than
% one cycle in the window is a trend (its period is not observed), which would let
% a sinusoid fit a trending series as if it were periodic (and below 1/(2N) its
% amplitude and phase are not determined: it cannot be told apart from a constant
% plus a linear trend, and the best-fitting frequency of a trending series tends to
% zero). The search is
% deterministic (no random starts: a grid search refined by Gauss-Newton iterations):
% see BF_FitSinusoids.
% The Fourier series are fitted by the Curve Fitting Toolbox from its own start
% points, which do not depend on the random number generator.

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
%% Fit the model
% ------------------------------------------------------------------------------
TSmodels = {'sin1', 'sin2', 'sin3', 'fourier1', 'fourier2', 'fourier3'}; % valid time-series models

if ~ismember(model, TSmodels)
	error('Invalid time-series model ''%s'' specified', model);
end

if size(y, 2) > size(y, 1)
	y = y';
end % y must be a column vector
N = length(y);

if strcmp(model(1:3), 'sin')
	% Sum of K sinusoids: least squares over amplitudes and phases, searched over frequencies
	K = str2double(model(4));
	numParams = 3*K; % amplitude, frequency and phase of each sinusoid
	if N <= numParams
		out = NaN; return
	end
	yfit = BF_FitSinusoids(y, K);
	res = y - yfit;
	sse = sum(res.^2);
	sstot = sum((y - mean(y)).^2);
	dfe = N - numParams; % degrees of freedom of the error
	out.r2 = 1 - sse/sstot;
	out.adjr2 = 1 - (1 - out.r2)*(N - 1)/dfe;
	out.rmse = sqrt(sse/dfe);
else
	% Fourier series (Curve Fitting Toolbox; not registered)
	BF_CheckToolbox('curve_fitting_toolbox');
	t = (1:N)'; % Time variable for equal sampling of the univariate time series
	try
		[~, gof, output] = fit(t, y, model); % fit the model
	catch emsg % this model can't even be fitted OR license problem
		if strcmp(emsg.message, 'NaN computed by model function.') || strcmp(emsg.message, 'Inf computed by model function.')
			fprintf(1, 'The model %s failed for this data -- returning NaNs for all fitting outputs\n', model);
			out = NaN; return
		else
			error('Unexpected error fitting ''%s'' to the time series', model)
		end
	end
	res = output.residuals;
	out.r2 = gof.rsquare;
	out.adjr2 = gof.adjrsquare;
	out.rmse = gof.rmse;
end

% ------------------------------------------------------------------------------
%% Remaining structure in the residuals
% ------------------------------------------------------------------------------
[out.resAC1, out.resAC2, out.resrunsz] = BF_ResidualStats(res, y);

end
